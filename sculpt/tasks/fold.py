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
    
    def __init__(self, model='Chai-1', num_structures=5, add_hydrogens=False, use_queue=False, scheduler=None, queue_args=None, scratch_dir=None):
        """Initialize the SculptFolder with model and number of structures.
        
        Args:
            model: The folding model to use ('Chai-1' or 'Boltz-2', default 'Chai-1').
            num_structures: The number of structures to generate per sequence (default 5).
            use_queue: Set to True to queue the jobs via a scheduler.
            scheduler: The name of the scheduler (e.g., 'SLURM', 'SGE') used if use_queue=True.
            queue_args: Additional arguments passed to queue() (e.g., node_name, queue, time).
                        Default: {'gpus': 1, 'mem': '16G'}
            scratch_dir: Custom directory to act as a persistent workspace instead of tempdir.
        """
        self.model = model
        self.add_hydrogens = add_hydrogens
        self.num_structures = num_structures
        self.use_queue = use_queue
        self.scheduler = scheduler
        self.queue_args = {'gpus': 1, 'mem': '16G'}
        if queue_args is not None:
            for key, value in queue_args.items():
                self.queue_args[key] = value
        self.scratch_dir = scratch_dir

    def fold(self, design: Design, ligand_sdf_files: List[str] = []):
        """Use Chai-1 or Boltz-2 to fold the FASTA file and optimize the ligand structure."""
        return self.fold_batch([design], ligand_sdf_files=ligand_sdf_files)[0]

    def fold_batch(self, designs: List[Design], ligand_sdf_files: List[str] = []):
        """Use Chai-1 or Boltz-2 to fold multiple FASTA files sequentially or concurrently.
        
        Args:
            designs: A list of Design objects to fold.
            ligand_sdf_files: A list of SDF files containing ligands to add to the structure.

        Returns:
            A list of lists, where each sublist contains the folded Design objects for a given input design.
        """
        import contextlib
        import os

        # We only accept 1 ligand SDF for now:
        if ligand_sdf_files is not None and len(ligand_sdf_files) > 0:
            ligand_sdf = ligand_sdf_files[0]
            # Convert to SMILES:
            ligand = StandardMolecule(structure_file=ligand_sdf)
            ligand_smiles = ligand.get_smiles()
        else:
            ligand_smiles = None

        if self.scratch_dir is not None:
            Path(self.scratch_dir).mkdir(parents=True, exist_ok=True)
            dir_context = contextlib.nullcontext(self.scratch_dir)
        else:
            dir_context = tempfile.TemporaryDirectory()

        with dir_context as work_dir:
            work_dir = Path(work_dir)
            job_ids = []
            tasks = []
            output_dirs = []

            for design in designs:
                design_dir = work_dir / design.name
                design_dir.mkdir(parents=True, exist_ok=True)
                
                fasta_file = design_dir / f"{design.name}.fasta"
                design.sequence_file(fasta_file)
                
                output_dir = design_dir / "folded"
                output_dirs.append(output_dir)

                if self.model == 'Chai-1':
                    task = ribbon.Chai1(
                            fasta_file=str(fasta_file),
                            output_dir=str(output_dir),
                            smiles_string=ligand_smiles,
                            num_ligands=1,
                            device='gpu'
                        )
                elif self.model == 'Boltz-2':
                    smiles_list = [ligand_smiles] if ligand_smiles else []
                    task = ribbon.Boltz2(
                        fasta_file=str(fasta_file),
                        output_dir=str(output_dir),
                        smiles_list=smiles_list,
                        device='gpu'
                    )
                else:
                    raise ValueError(f"Unknown folding model: {self.model}")
                
                tasks.append(task)
                if self.use_queue:
                    job_ids.append(task.queue(scheduler=self.scheduler, **self.queue_args))
                else:
                    task.run()

            # Wait for all folding jobs to finish
            if self.use_queue and job_ids:
                ribbon.wait_for_jobs(job_ids, scheduler=self.scheduler, max_wait=60*60*24)

            # Gather output files
            all_output_files = []
            for design, output_dir in list(zip(designs, output_dirs)):
                if self.model == 'Chai-1':
                    output_files = sorted(output_dir.glob('*_idx_*.cif'))
                elif self.model == 'Boltz-2':
                    output_files = sorted(output_dir.rglob('*_model_*.cif'))
                    if not output_files:
                        output_files = sorted(output_dir.rglob('*_model_*.pdb'))
                
                if len(output_files) > self.num_structures:
                    output_files = output_files[:self.num_structures]
                    
                all_output_files.append(output_files)
                print(f"[{design.name}] Collected {len(output_files)} structures.")

            # Process add_hydrogens
            if self.add_hydrogens:
                reduce_job_ids = []
                reduce_tasks = []
                reduce_out_files = []
                flattened_indices = []

                # Setup Reduce tasks
                for i, output_files in enumerate(all_output_files):
                    for j, file in enumerate(output_files):
                        idx = f"{i}_{j}"
                        task, target_outfile = self._create_reduce_task(
                            structure_file=file,
                            work_dir=work_dir / designs[i].name,
                            file_idx=idx,
                            ligand_sdf_files=ligand_sdf_files
                        )
                        flattened_indices.append((i, j, target_outfile))
                        if self.use_queue:
                            reduce_job_ids.append(task.queue(scheduler=self.scheduler, **self.queue_args))
                        else:
                            task.run()
                
                # Wait for Reduce jobs
                if self.use_queue and reduce_job_ids:
                    ribbon.wait_for_jobs(reduce_job_ids, scheduler=self.scheduler)
                    
                # Replace with hydrogenated structures
                for i, j, target_outfile in flattened_indices:
                    all_output_files[i][j] = target_outfile

            # Finalize designs
            all_folded_designs = []
            for design, output_files in list(zip(designs, all_output_files)):
                folded_designs = []
                for i, file in enumerate(output_files):
                    output_name = f"{design.name}_fold_{i}"
                    fold_design = design.copy(set_parent=True)
                    fold_design.name = output_name
                    fold_design.load_structure(file)
                    fold_design.history.append(f"Folded using {self.model}.")
                    folded_designs.append(fold_design)
                all_folded_designs.append(folded_designs)

        return all_folded_designs

    def _create_reduce_task(self, structure_file, work_dir, file_idx, ligand_sdf_files=[]):
        """Helper to create a Reduce Ribbon task and standard inputs/outputs within a specific directory."""
        infile = work_dir / f"reduce_in_{file_idx}.pdb"
        outfile = work_dir / f"reduce_out_{file_idx}.pdb"
        
        design = Design(name="temp")
        design.load_structure(structure_file)
        design.structure_file(str(infile))
        
        ligand_mols_sdfs = [(StandardMolecule(structure_file=sdf), sdf) for sdf in ligand_sdf_files]
        ligand_resnames_sdfs = []
        structure_hetero_residues = design.structure.get_residues(hetero_only=True)
        
        for mol, sdf in ligand_mols_sdfs:
            matched_residues = design.structure.match_residues_to_mol(structure_hetero_residues, mol)
            for match in matched_residues:
                chain, resnum, resname = match
                ligand_resnames_sdfs.append((resname, sdf))
        
        task = ribbon.Reduce(
            pdb_input_file=str(infile),
            pdb_output_file=str(outfile),
            custom_ligands=ligand_resnames_sdfs
        )
        return task, str(outfile)

    def add_hydrogens_to_file(self, structure_file, ligand_sdf_files=[]):
        """Add hydrogens to a structure file using Reduce inline. (Legacy method)"""
        with tempfile.TemporaryDirectory() as temp_dir:
            task, outfile = self._create_reduce_task(
                structure_file=structure_file,
                work_dir=Path(temp_dir),
                file_idx="legacy",
                ligand_sdf_files=ligand_sdf_files
            )
            task.run()
            # Need to copy to a safe place since temp_dir will be blown away
            out_fd, safe_outfile = tempfile.mkstemp(suffix='.pdb')
            os.close(out_fd)
            import shutil
            shutil.copy(outfile, safe_outfile)
            return safe_outfile

    def __call__(self, design: Design, ligand_sdf_files: List[str] = []):
        """Allow the Folder to be called directly.

        Args:
            design: The Design object to fold.
            ligand_sdf_files: Optional list of SDF files containing ligands to add hydrogens to.
        """
        return self.fold(design, ligand_sdf_files=ligand_sdf_files)


class DummyFolder:
    """Dummy Folder class for testing, substituting Chai-1 generation with fetching an existing mock design."""
    
    def __init__(self, model='Chai-1', num_structures=5, add_hydrogens=False, data_dir=None):
        """Initialize the DummyFolder.
        
        Args:
            model: The folding model to use (default 'Chai-1').
            num_structures: The number of structures to generate per sequence (default 5).
        """
        self.model = model
        self.add_hydrogens = add_hydrogens
        self.num_structures = num_structures
        
        # Path to the specific mock data folder
        self.data_dir = Path(data_dir)

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

    def fold_batch(self, designs: List[Design], ligand_sdf_files: List[str] = []):
        """Allow the Folder to be called directly."""
        output_designs = []
        for design in designs:
            output_designs.append(self.fold(design, ligand_sdf_files=ligand_sdf_files))
        return output_designs

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
