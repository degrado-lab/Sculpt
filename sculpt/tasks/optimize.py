### The following are the optimizer classes, used to optimize a structure file using EasyMD.
from typing import List, Union

import tempfile
from pathlib import Path
import os
import ribbon

from sculpt.design import Design
from sculpt.id import generate_id
from shim import StandardMolecule

class SculptOptimizer:
    """Optimizer class for structure optimization using custom restraints.
    
    This class provides functionality to optimize molecular structures using EasyMD
    with custom bonds, angles, and torsions as restraints.
    """
    
    def __init__(self, custom_bonds=[], custom_angles=[], custom_torsions=[], full_sim: bool = False, use_queue=False, scheduler=None, queue_args=None, scratch_dir=None):
        """Initialize the SculptOptimizer with custom restraints.
        
        Args:
            custom_bonds: List of custom bond restraints to apply during optimization.
            custom_angles: List of custom angle restraints to apply during optimization.
            custom_torsions: List of custom torsion restraints to apply during optimization.
            full_sim: Not implemented. If True, run a full simulation instead of just minimization (1 ns). This can be useful to "shakeout" the structure to better satisfy difficult restraints.
            use_queue: Set to True to queue the jobs via a scheduler.
            scheduler: The name of the scheduler (e.g., 'SLURM', 'SGE') used if use_queue=True.
            queue_args: Additional arguments passed to queue() (e.g., node_name, queue, time).
            scratch_dir: Custom directory to act as a persistent workspace instead of tempdir.
        """
        self.custom_bonds = custom_bonds
        self.custom_angles = custom_angles
        self.custom_torsions = custom_torsions
        self.full_sim = full_sim  # By default, just do minimization
        self.use_queue = use_queue
        self.scheduler = scheduler
        self.queue_args = {'gpus': 1, 'mem': '8G'}
        if queue_args is not None:
            for key, value in queue_args.items():
                self.queue_args[key] = value
        self.scratch_dir = scratch_dir

    def optimize(self, design: Design, sdf_files: List[str] = [], unique_naming: bool = False) -> Design:
        """Use EasyMD to optimize the structure file using custom restraints."""
        return self.optimize_batch([design], sdf_files=sdf_files, unique_naming=unique_naming)[0]
    def optimize_batch(self, designs: List[Design], sdf_files: List[str] = [], unique_naming: bool = False) -> List[Design]:
        """Use EasyMD to optimize a batch of structure files using custom restraints."""
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
            
            for design in designs:
                design_dir = work_dir / design.name
                design_dir.mkdir(parents=True, exist_ok=True)
                
                temp_structure_file = design_dir / f"{design.name}_opt_in.cif"
                design.structure_file(str(temp_structure_file))
                
                output_prefix = design_dir / "optimized"
                
                task = ribbon.EasyMD(
                    input_file=str(temp_structure_file),
                    output_prefix=str(output_prefix),
                    ligand_files=sdf_files,
                    duration=1,
                    custom_bonds=self.custom_bonds,
                    custom_torsions=self.custom_torsions,
                    custom_angles=self.custom_angles,
                    minimize_only= not self.full_sim
                )
                tasks.append(task)
                
                if self.use_queue:
                    job_ids.append(task.queue(scheduler=self.scheduler, **self.queue_args))
                else:
                    task.run()

            if self.use_queue and job_ids:
                ribbon.wait_for_jobs(job_ids, scheduler=self.scheduler)

            optimized_designs = []
            for design in designs:
                design_dir = work_dir / design.name
                
                ### DEBUG
                # print("Files in temp_dir:", os.listdir(design_dir))
                
                output_pdb_file = design_dir / "optimized_EM.pdb"

                if self.full_sim:
                    import mdtraj as md
                    output_pdb_file = design_dir / "optimized.pdb"
                    output_dcd_file = design_dir / "optimized_aligned.dcd"
                    if output_dcd_file.exists():
                        traj = md.load(str(output_dcd_file), top=str(output_pdb_file))
                        output_pdb_file = design_dir / "optimized_final_frame.pdb"
                        traj[-1].save_pdb(str(output_pdb_file))

                if unique_naming:
                    output_name = generate_id()
                else:
                    output_name = design.name + "_opt"

                optimized_design = design.copy(set_parent=True)
                optimized_design.name = output_name
                optimized_design.load_structure(output_pdb_file)

                standard_molecules = [StandardMolecule(structure_file=sdf) for sdf in sdf_files]
                optimized_design.structure.standardize(standard_molecules=standard_molecules, use_hydrogens=True, renumber=True)

                optimized_design.history.append("Optimized using EasyMD. \n"
                                                "Custom bonds: {}, \n"
                                                "Custom angles: {}, \n"
                                                "Custom torsions: {}".format(
                    self.custom_bonds, self.custom_angles, self.custom_torsions
                ))
                optimized_designs.append(optimized_design)

        return optimized_designs

    def __call__(self, structure_file, output_file, sdf_files = []):
        """Allow the optimizer to be called directly.
        
        Args:
            structure_file: The input structure file to optimize.
            output_file: The output file to write the optimized structure to.
            sdf_files: A list of SDF files to use for ligand parameterization.
            
        Returns:
            The result of the optimize method.
        """
        return self.optimize(structure_file, output_file, sdf_files) 



    
