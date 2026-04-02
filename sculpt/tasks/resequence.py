### The following are the optimizer classes, used to optimize a structure file using EasyMD.
import tempfile
from pathlib import Path
import ribbon

from sculpt.design import Design
from sculpt.id import generate_id

from shim.structure import ShimStructure

from typing import List

class SculptResequencer:
    """Resequencer class for generating new sequences from input structures.
    """
    
    def __init__(self, model='LigandMPNN', num_sequences=5, fixed_residues: str = None, use_queue=False, scheduler=None, queue_args=None, scratch_dir=None):
        """Initialize the SculptResequencer with model and number of sequences.
        
        Args:
            model: The folding model to use (default 'LigandMPNN').
            num_sequences: The number of sequences to generate (default 5).
            fixed_residues: A string specifying fixed residues (default None). (e.g. "A1 A2 A3")
            use_queue: Set to True to queue the jobs via a scheduler.
            scheduler: The name of the scheduler (e.g., 'SLURM', 'SGE') used if use_queue=True.
            queue_args: Additional arguments passed to queue() (e.g., node_name, queue, time).
            scratch_dir: Custom directory to act as a persistent workspace instead of tempdir.
        """
        if model not in ['LigandMPNN', 'LASErMPNN']:
            raise ValueError(f"Model '{model}' not recognized. Choose 'LigandMPNN' or 'LASErMPNN'.")
    
        self.model = model
        self.num_sequences = num_sequences
        self.fixed_residues = fixed_residues
        self.use_queue = use_queue
        self.scheduler = scheduler
        self.queue_args = queue_args or {}
        self.scratch_dir = scratch_dir

    def resequence(self, design: Design, unique_naming: bool = False):
        """Use LigandMPNN to generate new sequences from the input structure."""
        return self.resequence_batch([design], unique_naming=unique_naming)[0]

    def resequence_batch(self, designs: List[Design], unique_naming: bool = False):
        """Use LigandMPNN to generate new sequences from a batch of input structures."""
        import contextlib
        import os
        
        if self.scratch_dir is not None:
            Path(self.scratch_dir).mkdir(parents=True, exist_ok=True)
            dir_context = contextlib.nullcontext(self.scratch_dir)
        else:
            dir_context = tempfile.TemporaryDirectory()

        with dir_context as work_dir:
            work_dir = Path(work_dir)
            job_ids = []
            tasks = []

            if self.fixed_residues is not None:
                extra_args = f"--fixed_residues \"{self.fixed_residues}\""
            else:
                extra_args = ""

            for design in designs:
                design_dir = work_dir / design.name
                design_dir.mkdir(parents=True, exist_ok=True)
                input_structure_file = design_dir / f"{design.name}_in.pdb"
                
                if self.model == 'LigandMPNN':
                    design.structure_file(str(input_structure_file))
                    task = ribbon.LigandMPNN(
                        structure_list=[str(input_structure_file)],
                        output_dir=str(design_dir),
                        num_designs=self.num_sequences,
                        extra_args=extra_args
                    )
                elif self.model == 'LASErMPNN':
                    all_residues = design.structure.get_residues()
                    for chain_id, res_num, _ in all_residues:
                        if self.fixed_residues is not None and f"{chain_id}{res_num}" in self.fixed_residues:
                            design.structure.set_b_factor(chain_id, res_num, 1.0)
                        else:
                            design.structure.set_b_factor(chain_id, res_num, 0.0)

                    design.structure_file(str(input_structure_file))
                    task = ribbon.LASErMPNN(
                        structure_list=[str(input_structure_file)],
                        output_dir=str(design_dir),
                        num_designs=self.num_sequences,
                        fix_beta= True if self.fixed_residues else False,
                    )
                
                tasks.append(task)
                if self.use_queue:
                    job_ids.append(task.queue(scheduler=self.scheduler, **self.queue_args))
                else:
                    task.run()

            if self.use_queue and job_ids:
                ribbon.wait_for_jobs(job_ids, scheduler=self.scheduler)
                
            all_resequenced_designs = []
            for design in designs:
                design_dir = work_dir / design.name
                
                if self.model == 'LASErMPNN':
                    # Convert to FASTAs using shim:
                    lasermpnn_output_dir = design_dir / f"{design.name}_in"
                    seqs_split_dir = design_dir / 'seqs_split'
                    seqs_split_dir.mkdir(exist_ok=True)

                    if lasermpnn_output_dir.exists():
                        for pdb_file in lasermpnn_output_dir.glob('*.pdb'):
                            new_structure = ShimStructure(structure_file=str(pdb_file))
                            fasta_file = seqs_split_dir / (pdb_file.stem + '.fasta')
                            new_structure.to_fasta(str(fasta_file))

                resequenced_designs = []
                output_fasta_temp_dir = design_dir / 'seqs_split'

                if output_fasta_temp_dir.exists():
                    for fasta_file in output_fasta_temp_dir.glob('*.fasta'):
                        if unique_naming:
                            output_name = generate_id()
                        else:
                            output_name = design.name + "_reseq"

                        reseq_design = design.copy(set_parent=True)
                        reseq_design.name = output_name
                        reseq_design.load_sequence(str(fasta_file))
                        reseq_design.history.append(f"Resequenced using {self.model}.")
                        resequenced_designs.append(reseq_design)
                
                all_resequenced_designs.append(resequenced_designs)

        return all_resequenced_designs

    def __call__(self, design: Design, unique_naming: bool = False):
        """Allow the Resequencer to be called directly.

        Args:
            design: The Design object to resequence.
        """
        return self.resequence(design, unique_naming=unique_naming)
