### The following are the optimizer classes, used to optimize a structure file using EasyMD.
import tempfile
from pathlib import Path
import ribbon

from typing import List
from sculpt.design import Design
import random

from shim import StandardMolecule

class SculptFolder:
    """Folder class for folding FASTA files and ligands using Chai-1 or Boltz-2.
    NB. Known bug: Ligands get different resnames and chain IDs after adding Hydrogens with Reduce.
    """
    
    def __init__(self, model='Chai-1', num_structures=5, add_hydrogens=False):
        """Initialize the SculptFolder with model and number of structures.
        
        Args:
            model: The folding model to use ('Chai-1' or 'Boltz-2', default 'Chai-1').
            num_structures: The number of structures to generate per sequence (default 5).
        """
        self.model = model
        self.add_hydrogens = add_hydrogens
        self.num_structures = num_structures

    def fold(self, design: Design, ligand_sdf_files: List[str] = []):
        """Use Chai-1 or Boltz-2 to fold the FASTA file and optimize the ligand structure.
        
        Args:
            fasta_file: The input FASTA file to fold.
            output_dir: The directory to write the folded structures to.
            ligand_smiles: Optional SMILES string of the ligand to include in folding.

        N.B. Is there an easy way to convert SDF to SMILES here? Would the trouble of using SMILES as input.
        """

        # We only accept 1 ligand SDF for now:
        if ligand_sdf_files is not None and len(ligand_sdf_files) > 0:
            ligand_sdf = ligand_sdf_files[0]
            # Convert to SMILES:
            ligand = StandardMolecule(structure_file = ligand_sdf)
            ligand_smiles = ligand.get_smiles()
        else:
            ligand_smiles = None

        # Create a temporary directory to hold the outputs
        with tempfile.TemporaryDirectory() as temp_dir, \
                design.temp_sequence_file() as fasta_file:

            if self.model == 'Chai-1':
                ribbon.Chai1(
                        fasta_file = fasta_file,   # A single input FASTA. If there are multiple sequences, they will be folded in the same structure.
                        output_dir = str(Path(temp_dir) / "folded"),               # Where the outputs will be stored
                        smiles_string = ligand_smiles,      # SMILES string of our ligand
                        num_ligands = 1,                    # How many copies of our ligand?
                        device = 'gpu'                      # Run on GPU (necessary for Chai-1)
                    ).run()
    
                # Grab the output structures (prefix "_model_X.pdb") and copy them to the output directory. Change prefix to "_0, _1," etc.
                output_files = sorted( (Path(temp_dir) / "folded").glob('*_idx_*.cif') )
                if len(output_files) > self.num_structures:
                    output_files = output_files[:self.num_structures]
            elif self.model == 'Boltz-2':
                smiles_list = [ligand_smiles] if ligand_smiles else []
                ribbon.Boltz2(
                    fasta_file = fasta_file,
                    output_dir = str(Path(temp_dir) / "folded"),
                    smiles_list = smiles_list,
                    device = 'gpu'
                ).run()
                
                output_files = sorted( (Path(temp_dir) / "folded").rglob('*_model_*.cif') )
                if not output_files:
                    output_files = sorted( (Path(temp_dir) / "folded").rglob('*_model_*.pdb') )
                if len(output_files) > self.num_structures:
                    output_files = output_files[:self.num_structures]
            else:
                raise ValueError(f"Unknown folding model: {self.model}")

            print(output_files)

            if self.add_hydrogens:
                output_files_no_h = output_files.copy()
                # Add hydrogens to each structure
                output_files = []
                for file in output_files_no_h:
                    file_with_h = self.add_hydrogens_to_file(file, ligand_sdf_files)
                    output_files.append(file_with_h)

            folded_designs = []
            for i, file in enumerate(output_files):
                # Create the output Design object:
                output_name = design.name + f"_fold_{i}"

                # Make the new Design object:
                fold_design = design.copy(set_parent=True)
                fold_design.name = output_name
                fold_design.load_structure(file)

                # Add to the history:
                fold_design.history.append(f"Folded using {self.model}.")
                # if self.add_hydrogens:
                #     fold_design.history.append("Hydrogens added using Reduce.")

                folded_designs.append(fold_design)

        
        return folded_designs

    def add_hydrogens_to_file(self, structure_file, ligand_sdf_files=[]):
        """Add hydrogens to a structure file using Reduce.

        Args:
            structure_file: The input structure file to add hydrogens to.
            sdf_files: Optional list of SDF files containing ligands to add hydrogens to.
        """
        # Create a temp file for the input
        with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as infile:
            # We have an input file, and use a Design object to quickly convert it to PDB:
            design = Design(name="temp")
            design.load_structure(structure_file)
            design.structure_file(infile.name)

            # output
            outfile = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')

            # Create standard molecules
            ligand_mols_sdfs = [(StandardMolecule(structure_file = sdf), sdf) for sdf in ligand_sdf_files]

            # What is the residue that corresponds to each ligand?
            ligand_resnames_sdfs = []
            structure_hetero_residues = design.structure.get_residues(hetero_only=True)
            print("Hetero residues in structure:", structure_hetero_residues)
            for mol, sdf in ligand_mols_sdfs:
                matched_residues = design.structure.match_residues_to_mol(structure_hetero_residues, mol)
                # matched residue is list of (chain, resnum, resname)
                for match in matched_residues:
                    chain, resnum, resname = match
                    ligand_resnames_sdfs.append( (resname, sdf) )
                    print(f"Matched ligand residue {resname} in chain {chain} resnum {resnum} to SDF {sdf}")
                if len(matched_residues) == 0:
                    print(f"Warning: No residues in structure matched to ligand SDF {sdf}. Hydrogens will not be added to this ligand.")

            # make list of ligand names (LG1, LG2, etc.) and their corresponding sdf files
            custom_ligands = ligand_resnames_sdfs
            
            ribbon.Reduce(
                pdb_input_file=infile.name,
                pdb_output_file=outfile.name,
                custom_ligands=custom_ligands
            ).run()
        
        return outfile.name

    def __call__(self, design: Design, ligand_sdf_files: List[str] = []):
        """Allow the Folder to be called directly.

        Args:
            design: The Design object to fold.
            ligand_sdf_files: Optional list of SDF files containing ligands to add hydrogens to.
        """
        return self.fold(design, ligand_sdf_files=ligand_sdf_files)


class DummyFolder:
    """Dummy Folder class for testing, substituting Chai-1 generation with fetching an existing mock design."""
    
    def __init__(self, model='Chai-1', num_structures=5, add_hydrogens=False):
        """Initialize the DummyFolder.
        
        Args:
            model: The folding model to use (default 'Chai-1').
            num_structures: The number of structures to generate per sequence (default 5).
        """
        self.model = model
        self.add_hydrogens = add_hydrogens
        self.num_structures = num_structures
        
        # Path to the specific mock data folder
        self.data_dir = Path(__file__).resolve().parent.parent.parent / "data" / "1OHP_mutants"

    def fold(self, design: Design, ligand_sdf_files: List[str] = []):
        """Pretend to fold the FASTA file by randomly picking existing CIF files.
        """
        
        print("Using DummyFolder to fold the FASTA file.")
        
        # We only accept 1 ligand SDF for now:
        if ligand_sdf_files is not None and len(ligand_sdf_files) > 0:
            ligand_sdf = ligand_sdf_files[0]
            # Convert to SMILES:
            ligand = StandardMolecule(structure_file = ligand_sdf)
            ligand_smiles = ligand.get_smiles()
        else:
            ligand_smiles = None
            
        available_files = list(self.data_dir.glob("*.cif"))
        if not available_files:
            raise FileNotFoundError(f"No CIF files found in {self.data_dir} for DummyFolder")
            
        # Select random structures for mock output
        output_files = random.choices(available_files, k=self.num_structures)
        
        folded_designs = []
        for i, file in enumerate(output_files):
            # Create the output Design object:
            output_name = design.name + f"_fold_{i}"

            # Make the new Design object:
            fold_design = design.copy(set_parent=True)
            fold_design.name = output_name
            fold_design.load_structure(str(file))

            # Add to the history:
            fold_design.history.append("Folded using DummyFolder mock.")

            folded_designs.append(fold_design)

        return folded_designs

    def __call__(self, design: Design, ligand_sdf_files: List[str] = []):
        """Allow the Folder to be called directly."""
        return self.fold(design, ligand_sdf_files=ligand_sdf_files)



class SculptHydrogenAdder:
    """Folder class for folding FASTA files and ligands using Chai-1.
    """
    
    def __init__(self):
        """Initialize the SculptFolder with model and number of structures.
        """
        pass

    def add_hydrogens(self, design: Design, ligand_sdf_files=[]):
        """Add hydrogens to a structure file using Reduce.

        Args:
            structure_file: The input structure file to add hydrogens to.
            sdf_files: Optional list of SDF files containing ligands to add hydrogens to.

        Returns:
            A new Design object with hydrogens added.
        """
        # Create a temp file for the input
        with tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as infile, \
             tempfile.NamedTemporaryFile(delete=False, suffix='.pdb') as outfile:
            
            # input
            design.structure_file(infile.name)

            # output
            outfile = tempfile.NamedTemporaryFile(delete=False, suffix='.pdb')

            # Create standard molecules
            ligand_mols_sdfs = [(StandardMolecule(structure_file = sdf), sdf) for sdf in ligand_sdf_files]

            # What is the residue that corresponds to each ligand?
            ligand_resnames_sdfs = []
            structure_hetero_residues = design.structure.get_residues(hetero_only=True)
            print("Hetero residues in structure:", structure_hetero_residues)
            for mol, sdf in ligand_mols_sdfs:
                matched_residues = design.structure.match_residues_to_mol(structure_hetero_residues, mol)
                # matched residue is list of (chain, resnum, resname)
                for match in matched_residues:
                    chain, resnum, resname = match
                    ligand_resnames_sdfs.append( (resname, sdf) )
                    print(f"Matched ligand residue {resname} in chain {chain} resnum {resnum} to SDF {sdf}")
                if len(matched_residues) == 0:
                    print(f"Warning: No residues in structure matched to ligand SDF {sdf}. Hydrogens will not be added to this ligand.")

            # make list of ligand names (LG1, LG2, etc.) and their corresponding sdf files
            custom_ligands = ligand_resnames_sdfs
            
            ribbon.Reduce(
                pdb_input_file=infile.name,
                pdb_output_file=outfile.name,
                custom_ligands=custom_ligands
            ).run()
            
            # Make the new Design object:
            h_design = design.copy(set_parent=True)
            h_design.name = design.name
            h_design.load_structure(outfile.name)

            # Add to the history:
            h_design.history.append("Hydrogens added using Reduce.")

        return h_design

    def __call__(self, design: Design, ligand_sdf_files: List[str] = []):
        """Allow the Folder to be called directly.

        Args:
            design: The Design object to fold.
            ligand_sdf_files: Optional list of SDF files containing ligands to add hydrogens to.
        """
        return self.add_hydrogens(design, ligand_sdf_files=ligand_sdf_files)
