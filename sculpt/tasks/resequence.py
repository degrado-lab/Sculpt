### The following are the optimizer classes, used to optimize a structure file using EasyMD.
import tempfile
from pathlib import Path
import ribbon

from sculpt.design import Design
from sculpt.id import generate_id

from shim.structure import ShimStructure

class SculptResequencer:
    """Resequencer class for generating new sequences from input structures.
    """
    
    def __init__(self, model='LigandMPNN', num_sequences=5, fixed_residues: str = None):
        """Initialize the SculptResequencer with model and number of sequences.
        
        Args:
            model: The folding model to use (default 'LigandMPNN').
            num_sequences: The number of sequences to generate (default 5).
            fixed_residues: A string specifying fixed residues (default None). (e.g. "A1 A2 A3")
        """
        if model not in ['LigandMPNN', 'LASErMPNN']:
            raise ValueError(f"Model '{model}' not recognized. Choose 'LigandMPNN' or 'LASErMPNN'.")
    
        self.model = model
        self.num_sequences = num_sequences
        self.fixed_residues = fixed_residues

    def resequence(self, design: Design, unique_naming: bool = False):
        """Use LigandMPNN to generate new sequences from the input structure.

        Args:
            design: The Design object to generate sequences from.
        """
        with tempfile.TemporaryDirectory() as temp_dir:

            if self.fixed_residues is not None:
                extra_args = f"--fixed_residues \"{self.fixed_residues}\""
            else:
                extra_args = ""

            if self.model == 'LigandMPNN':
                with design.temp_structure_file() as input_structure_file:
                    print(f"Running LigandMPNN on {input_structure_file}...")
                    ribbon.LigandMPNN(
                        structure_list=[str(input_structure_file)],
                        output_dir=temp_dir,
                        num_designs=self.num_sequences,
                        extra_args=extra_args
                    ).run()
            
            elif self.model == 'LASErMPNN':
                # First, we create 5 new sequences for this structure:

                # Make a copy of the structure, where we set the b-factors to 1.0 for all fixed residues:
                all_residues = design.structure.get_residues()
                for chain_id, res_num, _ in all_residues:
                    if self.fixed_residues is not None and f"{chain_id}{res_num}" in self.fixed_residues:
                        design.structure.set_b_factor(chain_id, res_num, 1.0)
                    else:
                        design.structure.set_b_factor(chain_id, res_num, 0.0)

                with design.temp_structure_file(suffix='.pdb') as input_structure_file:
                    print(f"Running LASErMPNN on {input_structure_file}...")
                    ribbon.LASErMPNN(
                        structure_list=[str(input_structure_file)],
                        output_dir=temp_dir,
                        num_designs=self.num_sequences,
                        fix_beta= True if self.fixed_residues else False,
                    ).run()
                
                    # This will output PDBs. Convert to FASTAs using shim:
                    lasermpnn_output_dir = Path(temp_dir) / input_structure_file.stem

                    # list files in directory:
                    for file in lasermpnn_output_dir.iterdir():
                        print(' -', file.name)

                    # Make seqs_split directory:
                    (Path(temp_dir) / 'seqs_split').mkdir(exist_ok=True)

                    for pdb_file in lasermpnn_output_dir.glob('*.pdb'):
                        new_structure = ShimStructure(structure_file=pdb_file)
                        fasta_file = Path(temp_dir) / 'seqs_split' / (pdb_file.stem + '.fasta')
                        print(f"Converting {pdb_file} to {fasta_file}...")
                        new_structure.to_fasta(fasta_file)

            # Sequences are in the [tempfile]/sequences_split/ directory. 
            # We'll create a new Design for each sequence:
            resequenced_designs = []
            output_fasta_temp_dir = Path(temp_dir) / 'seqs_split'

            for fasta_file in output_fasta_temp_dir.glob('*.fasta'):
                # Create the output Design object:
                if unique_naming:
                    output_name = generate_id()
                else:
                    output_name = design.name + "_reseq"

                # Make the new Design object:
                reseq_design = design.copy(set_parent=True)
                reseq_design.name = output_name
                reseq_design.load_sequence(fasta_file)

                # Add to the history:
                reseq_design.history.append("Resequenced using LigandMPNN.")

                resequenced_designs.append(reseq_design)

        return resequenced_designs

    def __call__(self, design: Design, unique_naming: bool = False):
        """Allow the Resequencer to be called directly.

        Args:
            design: The Design object to resequence.
        """
        return self.resequence(design, unique_naming=unique_naming)
