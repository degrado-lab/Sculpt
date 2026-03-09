"""
Crossover strategies for genetic algorithms using Design objects.

This file defines a structural Protocol `CrossoverStrategy` and a concrete
implementation: SimpleCrossover.
"""

from typing import Sequence, Protocol, List, TypeVar, Tuple
import random
import numpy as np
from Bio.PDB import Superimposer, MMCIFParser, PDBParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from sculpt.id import generate_id

Design = TypeVar('Design')


class CrossoverStrategy(Protocol[Design]):
    """Abstract crossover strategy for genetic algorithms.

    Implementations should return a sequence of offspring individuals from the
    provided structured population.
    """

    def crossover(self, population: Sequence[Design], **kwargs) -> List[Design]:
        pass


#TODO: Move helper functions to utils submodule?
def _parse_fasta(fasta_text: str) -> Tuple[str, str]:
    """
    Helper to parse a simple FASTA string into header and sequence.
    Assumes a single object with a single header."""
    if not fasta_text:
        return "", ""
    lines = fasta_text.strip().split('\n')
    if not lines:
        return "", ""
    header = lines[0] if lines[0].startswith('>') else ">"
    # The sequence could be on multiple lines or just one
    seq = "".join([line.strip() for line in lines if not line.startswith('>')])
    if not seq and not lines[0].startswith('>'):
        seq = lines[0].strip()
    return header, seq


def _create_fasta(header: str, seq: str, line_length: int = 80) -> str:
    """Helper to create a FASTA string."""
    header = header if header.startswith('>') else f">{header}"
    seq_lines = [seq[i:i+line_length] for i in range(0, len(seq), line_length)]
    return f"{header}\n" + "\n".join(seq_lines)


class SimpleCrossover:
    """Simple one-point sequence crossover for a population of Design objects.

    This crossover strategy randomly pairs individuals from the population and
    performs a single-point crossover on their sequences with a given probability.
    Structure data (if any) is not crossed over; only the sequence is modified.

    Args:
        crossover_rate: Probability of crossover occurring for a pair.
    """

    def __init__(self, crossover_rate: float = 0.5):
        self.crossover_rate = max(0.0, min(1.0, float(crossover_rate)))

    def crossover(self, population: Sequence[Design], unique_naming: bool = True) -> List[Design]:
        """Perform crossover on a population.
        
        Args:
            population: Sequence of Design objects to cross over.
            unique_naming: If True, generate unique names for the offspring.
                           Otherwise, derive names from parents.
        
        Returns:
            A new list of offspring Design objects.
        """
        n = len(population)
        if n == 0:
            return []
        if n == 1:
            return [population[0].copy()]

        # Shuffle the population to create random pairs
        indices = list(range(n))
        random.shuffle(indices)
        
        offspring_pop = []
        
        # Process in pairs
        for i in range(0, n, 2):
            parent1 = population[indices[i]]
            
            if i + 1 < n:
                parent2 = population[indices[i+1]]
            else:
                # Odd number of individuals, just copy the last one without crossover
                offspring_pop.append(parent1.copy())
                continue
                
            # Create offspring
            child1 = parent1.copy(set_parent=True)
            child1.parents = [parent1, parent2]
            
            child2 = parent2.copy(set_parent=True)
            child2.parents = [parent2, parent1]
            
            # Apply crossover based on rate
            if random.random() < self.crossover_rate and child1.sequence and child2.sequence:
                h1, seq1 = _parse_fasta(child1.sequence)
                h2, seq2 = _parse_fasta(child2.sequence)
                
                # We can only perform sequence crossover if both sequences have length >= 1
                min_len = min(len(seq1), len(seq2))
                
                if min_len > 1:
                    # Pick a random crossover point (between 1 and min_len - 1)
                    cx_point = random.randint(1, min_len - 1)
                    
                    new_seq1 = seq1[:cx_point] + seq2[cx_point:]
                    new_seq2 = seq2[:cx_point] + seq1[cx_point:]
                    
                    if unique_naming:
                        c1_name = generate_id()
                        c2_name = generate_id()
                    else:
                        c1_name = f"{parent1.name}_x_{parent2.name}"
                        c2_name = f"{parent2.name}_x_{parent1.name}"
                        
                    child1.name = c1_name
                    child2.name = c2_name
                    
                    # Create new fasta text
                    child1.sequence = _create_fasta(c1_name, new_seq1)
                    child2.sequence = _create_fasta(c2_name, new_seq2)

                    # Because we modified the sequence, the structure is no longer strictly matching it.
                    # We should probably clear it since fold step will regenerate it.
                    child1.structure = None
                    child2.structure = None
                    
                    child1.history.append(f"Crossover with {parent2.name} at position {cx_point}")
                    child2.history.append(f"Crossover with {parent1.name} at position {cx_point}")
                    
            offspring_pop.extend([child1, child2])
            
        return offspring_pop

    def __call__(self, population: Sequence[Design], unique_naming: bool = False) -> List[Design]:
        """Allow the crossover strategy to be called directly.
        
        Args:
            population: Sequence of Design objects to cross over.
            unique_naming: If True, generate unique names for the offspring.
            
        Returns:
            The result of the crossover method.
        """
        return self.crossover(population, unique_naming=unique_naming)



class SpatialCrossover:
    """ 3D crossover algorithm for a population of Design objects.

    This crossover strategy randomly pairs individuals from the population and
    aligns their structures. 
    Two residues (A and B) are chosen from parent 1. The first child receives 
    the residues from parent 1 which are closest to point A, and the residues 
    from parent 2 which are closer to B. The second child receives the opposite.
    
    Args:
        crossover_rate: Probability of crossover occurring for a pair.
    """

    def __init__(self, crossover_rate: float = 0.5):
        self.crossover_rate = max(0.0, min(1.0, float(crossover_rate)))

    def crossover(self, population: Sequence[Design], unique_naming: bool = True) -> List[Design]:
        """Perform crossover on a population.
        
        Args:
            population: Sequence of Design objects to cross over.
            unique_naming: If True, generate unique names for the offspring.
                           Otherwise, derive names from parents.
        
        Returns:
            A new list of offspring Design objects.
        """
        n = len(population)
        if n == 0:
            return []
        if n == 1:
            return [population[0].copy()]

        # Shuffle the population to create random pairs
        indices = list(range(n))
        random.shuffle(indices)
        
        offspring_pop = []
        
        # Process in pairs
        for i in range(0, n, 2):
            parent1 = population[indices[i]]
            
            if i + 1 < n:
                parent2 = population[indices[i+1]]
            else:
                # Odd number of individuals, just copy the last one without crossover
                offspring_pop.append(parent1.copy())
                continue
                
            # Apply crossover based on rate and availability of structures
            if random.random() < self.crossover_rate and parent1.structure and parent2.structure:
                try:
                    # 1. Parse structures using Bio.PDB
                    mcif_parser = MMCIFParser(QUIET=True)
                    pdb_parser = PDBParser(QUIET=True)
                    
                    def parse_struct(design, name):
                        with design.temp_structure_file() as f:
                            if f.suffix.lower() in ['.cif', '.mmcif']:
                                return mcif_parser.get_structure(name, str(f))
                            else:
                                return pdb_parser.get_structure(name, str(f))

                    struct1 = parse_struct(parent1, "p1")
                    struct2 = parse_struct(parent2, "p2")
                    
                    # 2. Align struct2 to struct1
                    self._align_by_ca(struct1, struct2)
                    
                    # 3. Define points A and B (e.g., N and C termini of parent 1)
                    point_a = self._get_res_ca_coords(struct1, residue_id=1)
                    point_b = self._get_res_ca_coords(struct1, residue_id=-1)
                    
                    # 4. Determine residency identity for each residue based on proximity to A or B
                    identity = self._get_residue_a_b_identity(struct1, point_a, point_b)

                    # 5. Create children sequences
                    res1 = [r for r in struct1.get_residues() if r.id[0] == " "]
                    res2 = [r for r in struct2.get_residues() if r.id[0] == " "]
                    
                    min_res = min(len(res1), len(res2), len(identity))
                    
                    residues_c1 = []
                    residues_c2 = []
                    for idx in range(min_res):
                        if identity[idx] == "A":
                            residues_c1.append(res1[idx])
                            residues_c2.append(res2[idx])
                        else:
                            residues_c1.append(res2[idx])
                            residues_c2.append(res1[idx])
                    
                    # Create offspring objects
                    child1 = parent1.copy(set_parent=True)
                    child2 = parent2.copy(set_parent=True)
                    child1.parents = [parent1, parent2]
                    child2.parents = [parent2, parent1]
                    
                    if unique_naming:
                        c1_name = generate_id()
                        c2_name = generate_id()
                    else:
                        c1_name = f"{parent1.name}_spatial_{parent2.name}"
                        c2_name = f"{parent2.name}_spatial_{parent1.name}"
                        
                    child1.name = c1_name
                    child2.name = c2_name
                    
                    # Create new fasta text
                    seq1_str = self._residues_to_seq(residues_c1)
                    seq2_str = self._residues_to_seq(residues_c2)
                    
                    child1.sequence = _create_fasta(c1_name, seq1_str)
                    child2.sequence = _create_fasta(c2_name, seq2_str)

                    # Clear structures as they need folding for the new sequences
                    child1.structure = None
                    child2.structure = None
                    
                    child1.history.append(f"Spatial crossover with {parent2.name}")
                    child2.history.append(f"Spatial crossover with {parent1.name}")
                    
                    offspring_pop.extend([child1, child2])
                    continue
                except Exception as e:
                    # In case of any error in spatial processing, fall back to simple copies
                    pass
            
            # Default fallback: Copy parents if crossover not performed or failed
            offspring_pop.extend([parent1.copy(), parent2.copy()])
            
        return offspring_pop

    def __call__(self, population: Sequence[Design], unique_naming: bool = False) -> List[Design]:
        """Allow the crossover strategy to be called directly.
        
        Args:
            population: Sequence of Design objects to cross over.
            unique_naming: If True, generate unique names for the offspring.
            
        Returns:
            The result of the crossover method.
        """
        return self.crossover(population, unique_naming=unique_naming)

    def _get_res_ca_coords(self, structure, model_id=0, chain_id="A", residue_id=1):
        """Helper to get CA coordinates of a specific residue."""
        model = structure[model_id]
        if chain_id not in model:
            # Fallback to the first chain if "A" is not found
            chain = list(model.get_chains())[0]
        else:
            chain = model[chain_id]
            
        residues = [r for r in chain if r.id[0] == " "]
        if not residues:
            raise ValueError("No standard residues found in chain")
            
        if residue_id < 0:
            target_res = residues[residue_id]
        else:
            target_res = None
            for res in residues:
                if res.id[1] == residue_id:
                    target_res = res
                    break
            if target_res is None:
                target_res = residues[0]
                
        if "CA" not in target_res:
            raise ValueError(f"No CA atom in residue {target_res.id}")
        return target_res["CA"].coord

    def _get_ca_atoms(self, structure, model_id=0, chain_id="A"):
        """Return dict {(resseq, icode): CA_atom} for a chain."""
        model = structure[model_id]
        if chain_id not in model:
            chain = list(model.get_chains())[0]
        else:
            chain = model[chain_id]
            
        ca = {}
        for res in chain:
            if res.id[0] != " " or "CA" not in res:
                continue
            ca[(res.id[1], res.id[2])] = res["CA"]
        return ca

    def _align_by_ca(self, ref, mob, chain_ref="A", chain_mobile="A", model_id=0):
        """Align mobile structure to reference structure using CA atoms."""
        ref_ca = self._get_ca_atoms(ref, model_id=model_id, chain_id=chain_ref)
        mob_ca = self._get_ca_atoms(mob, model_id=model_id, chain_id=chain_mobile)
        
        common_keys = sorted(set(ref_ca.keys()) & set(mob_ca.keys()))
        if len(common_keys) < 3:
            raise ValueError(f"Need >=3 matched CA atoms; found {len(common_keys)}")
            
        fixed_atoms = [ref_ca[k] for k in common_keys]
        moving_atoms = [mob_ca[k] for k in common_keys]
        
        sup = Superimposer()
        sup.set_atoms(fixed_atoms, moving_atoms)
        sup.apply(mob.get_atoms())
        return sup.rms

    def _get_residue_a_b_identity(self, structure_a, point_a, point_b):
        """Determine for each residue if it's closer to point A or point B."""
        identity = []
        for res in structure_a.get_residues():
            if res.id[0] != " ":
                continue
            if "CA" in res:
                ca = res["CA"]
                dist_a = np.linalg.norm(ca.coord - point_a)
                dist_b = np.linalg.norm(ca.coord - point_b)
                if dist_a < dist_b:
                    identity.append("A")
                else:
                    identity.append("B")
            else:
                identity.append("A") # Default fallback
        return identity

    def _residues_to_seq(self, residues):
        """Convert a list of Bio.PDB residues to a sequence string."""
        three_to_one = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        return "".join([three_to_one.get(res.resname, 'X') for res in residues])