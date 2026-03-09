from Bio.PDB import MMCIFParser, Superimposer, PDBIO, Select
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

### V1 Spatial Crossover Algorithm, Nicholas Freitas 2/6/2026
### Super quick and dirty. First, takes in two structures and aligns them. User provides 2 points
### in space, and for each residue in structure A, checks which point it is closer to. Swaps
### the sequence from structure B to create a new sequence, output as fasta.


### Helper functions
def get_res_ca_coords(structure, model_id=0, chain_id="A", residue_id=1):
    model = structure[model_id]
    chain = model[chain_id]

    # if the residue_id is negative, count from the end
    if residue_id < 0:
        residue_id = len(chain) + residue_id + 1 # resid is typically 1-indexed

    for res in chain:
        if res.id[0] != " ":  # skip HETATM/water
            continue
        if res.id[1] != residue_id:
            continue
        return res["CA"].coord
    raise ValueError("Residue not found")

def _iter_standard_residues(structure, model_id=0, chain_id="A"):
    """Yield standard (non-hetero) residues from a given chain."""
    model = structure[model_id]
    chain = model[chain_id]
    for res in chain:
        if res.id[0] != " ":  # skip HETATM/water
            continue
        yield res

def _residue_atom_positions(residue, heavy_only=True):
    """Return (N,3) positions for atoms in residue."""
    pts = []
    for atom in residue.get_atoms():
        if heavy_only and atom.element == "H":
            continue
        pts.append(atom.coord)
    if not pts:
        return np.zeros((0, 3), dtype=float)
    return np.asarray(pts, dtype=float)

def _ca_centroid(structure, model_id=0, chain_id="A"):
    """Centroid of CA atoms for standard residues (not mass-weighted)."""
    coords = []
    for res in _iter_standard_residues(structure, model_id, chain_id):
        if "CA" in res:
            coords.append(res["CA"].coord)
    if len(coords) == 0:
        raise ValueError("No CA atoms found for centroid calculation.")
    return np.mean(np.asarray(coords, dtype=float), axis=0)


### Alignment

def get_ca_atoms(structure, model_id=0, chain_id="A"):
    """Return dict {(resseq, icode): CA_atom} for a chain."""
    model = structure[model_id]
    chain = model[chain_id]
    ca = {}
    for res in chain:
        # Skip hetero/waters; keep standard residues
        if res.id[0] != " ":
            continue
        if "CA" not in res:
            continue
        resseq, icode = res.id[1], res.id[2]
        ca[(resseq, icode)] = res["CA"]
    return ca

def align_by_ca(ref, mob, output_file=None,
                chain_ref="A", chain_mobile="A", model_id=0):

    ref_ca = get_ca_atoms(ref, model_id=model_id, chain_id=chain_ref)
    mob_ca = get_ca_atoms(mob, model_id=model_id, chain_id=chain_mobile)

    # Match residues by (resseq, insertion_code)
    common_keys = sorted(set(ref_ca.keys()) & set(mob_ca.keys()))
    if len(common_keys) < 3:
        raise ValueError(f"Need >=3 matched CA atoms; found {len(common_keys)}")

    fixed_atoms = [ref_ca[k] for k in common_keys]
    moving_atoms = [mob_ca[k] for k in common_keys]

    sup = Superimposer()
    sup.set_atoms(fixed_atoms, moving_atoms)  # computes rotation+translation

    # Apply transform to *all atoms* in the mobile structure
    sup.apply(mob.get_atoms())

    print(f"RMSD over {len(common_keys)} CA atoms: {sup.rms:.3f} Å")

    # Write aligned mobile structure
    if output_file:
        io = PDBIO()
        io.set_structure(mob)
        io.save(output_file)

    return ref, mob

### Spatial Crossover

def get_residue_a_b_identity(structure_a, point_a, point_b):
    '''For each residue in structure A, check if it's closer to point A or point B
    '''

    identity = []
    for res in structure_a.get_residues():
        ca = res["CA"]
        if np.linalg.norm(ca.coord - point_a) < np.linalg.norm(ca.coord - point_b):
            identity.append("A")
        else:
            identity.append("B")
    return identity
    
def get_crossover_sequence(structure_a, structure_b, identity):
    '''Given two structures and their identity (A or B), return the crossover sequence
    '''

    crossover_sequence = []
    for i in range(len(identity)):
        if identity[i] == "A":
            crossover_sequence.append(structure_a[0]["A"][(" ", i + 1, " ")])
        else:
            crossover_sequence.append(structure_b[0]["A"][(" ", i + 1, " ")])

    # Create fasta:



    return crossover_sequence
    


def save_residues_to_fasta(residues, output_file="output.fasta"):
    """Converts a list of Bio.PDB residues to a FASTA file."""
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    sequence = ""
    for res in residues:
        sequence += three_to_one.get(res.resname, 'X')

    record = SeqRecord(
        Seq(sequence),
        id="crossover_chimera",
        description=""
    )

    with open(output_file, "w") as handle:
        SeqIO.write(record, handle, "fasta")
    print(f"FASTA saved to {output_file}")

##################################################
####### Running the algorithm ####################
##################################################

structure_file_a = "HG3_idx_0.cif"
structure_file_b = "HG3_idx_1.cif"

parser = MMCIFParser(QUIET=True)
ref = parser.get_structure("ref", structure_file_a)
mob = parser.get_structure("mob", structure_file_b)

# Load the files, and return the aligned structures:
struct_aligned_a, struct_aligned_b = align_by_ca(ref, mob, chain_ref="A", chain_mobile="A")

# Establish two points in space (nM)
# Here, I'm choosing the N and C terminus, but this is arbitrary.
point_a = get_res_ca_coords(struct_aligned_a, residue_id=1)
point_b = get_res_ca_coords(struct_aligned_b, residue_id=-1)

# For each residue in structure A, check if it's closer to point A or point B (it's "parent")
# This is just referencing structure A
residue_parent_list = get_residue_a_b_identity(struct_aligned_a, point_a, point_b)

# Get the crossover sequence
crossover_sequence = get_crossover_sequence(struct_aligned_a, struct_aligned_b, residue_parent_list)

print(crossover_sequence)
save_residues_to_fasta(crossover_sequence, "crossover.fasta")