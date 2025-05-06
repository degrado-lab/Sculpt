import json
import shutil
from pathlib import Path
import ribbon
from ribbon.utils import make_directories, list_files


def copy_input_files(input_dir: Path, target_dir: Path):
    """
    Copy all files from the input directory to the target directory.
    """
    for file in input_dir.iterdir():
        shutil.copy(file, target_dir)

def energy_minimization(previous_dir: Path, em_dir: Path, sdf_file: str, custom_bonds=[], custom_torsions=[], custom_angles=[], scheduler: str = 'local', queue: str = '*', node_name: str = '*'):
    """
    For every structure in the previous cycle, calculate the distance between
    ligand and protein atoms then run energy minimization.
    """
    # Look for both .cif and .pdb files in the previous design directory.
    prev_design_dir = previous_dir / '5_Top_Designs'
    structure_files = list_files(prev_design_dir, '.cif') + list_files(prev_design_dir, '.pdb')
    print('Number of structures:', len(structure_files))
    
    # For each item in custom_bonds, etc, remove the whitespace.
    custom_bonds = [bond.replace(' ', '') for bond in custom_bonds]
    custom_torsions = [torsion.replace(' ', '') for torsion in custom_torsions]
    custom_angles = [angle.replace(' ', '') for angle in custom_angles]

    job_list = []
    for file in structure_files:

        output_prefix = em_dir / Path(file).stem

        # Run energy minimization with custom restraints
        EM_task = ribbon.EasyMD(
            input_file=file,
            output_prefix=output_prefix,
            ligand_files=[sdf_file],
            duration=1,
            custom_bonds=custom_bonds,
            custom_torsions=custom_torsions,
            custom_angles=custom_angles,
            minimize_only=True,)
        
        if scheduler == 'local':
            # Run locally
            EM_task.run()
        elif scheduler in ['SGE', 'SLURM']:
            # Queue the job using SGE
            job_id = EM_task.queue(
                queue=queue,
                output_file=str('./job_outputs') + f'/EM_{Path(file).stem}.out',
                scheduler=scheduler,
                time='00:30:00',
                job_name=f'EM_{Path(file).stem}',
                node_name = node_name) #'qb3-atgpu*'
            job_list.append(job_id)
            print('Queueing EM Job: ', job_id)
    
    # Wait for all jobs to finish
    if scheduler == 'local':
        # If running locally, we don't need to wait for jobs.
        print('All jobs completed.')
    elif scheduler in ['SGE', 'SLURM']:
        # Wait for all jobs to finish
        print('Waiting for all jobs to finish...')
        ribbon.wait_for_jobs(job_list, scheduler=scheduler, max_wait=60*60*24) # If we're not done in 24 hours, kill it.
        
def run_ligand_mpnn(em_dir: Path, LMPNN_dir: Path, design_num: int, residues_to_change: str, temperature: float = 0.1):
    """
    Run LigandMPNN to generate sequences for each energy-minimized structure.
    Also, fix the filenames in the 'seqs_split' folder.
    """
    # Gather EM-minimized structures that end with "_EM.pdb"
    EM_structures_list = list_files(em_dir, '_EM.pdb')
    extra_args = f'--redesigned_residues "{residues_to_change}" --homo_oligomer 1 --temperature {str(temperature)}'
    ribbon.LigandMPNN(
        LMPNN_dir,
        structure_list=EM_structures_list,
        num_designs=design_num,
        extra_args=extra_args
    ).run()

    # Rename the files to remove a possible '.cif' substring in the filename.
    seqs_split_dir = LMPNN_dir / 'seqs_split'
    for sequence_file in list_files(seqs_split_dir, '.fasta'):
        sequence_file = Path(sequence_file)
        new_name = sequence_file.name.replace('.cif', '_cif')
        shutil.move(sequence_file, sequence_file.parent / new_name)

def run_chai(LMPNN_dir: Path, Chai_dir: Path, smiles: str, num_ligands: int, scheduler: str = 'local', queue: str = '*', node_name: str = '*'):
    """
    For each designed sequence file, run CHAI-1 to predict structures.
    """
    seqs_split_dir = LMPNN_dir / 'seqs_split'
    chai_job_id_list = []
    for fasta_file in list_files(seqs_split_dir, '.fasta'):
        fasta_file = Path(fasta_file)
        Chai_task = ribbon.Chai1(
            fasta_file=fasta_file,
            smiles_string=smiles,
            num_ligands=num_ligands,
            output_dir=Chai_dir / fasta_file.stem,
        )
        
        if scheduler == 'local':
            # Run locally
            Chai_task.run()
        elif scheduler in ['SGE', 'SLURM']:
            # Queue the job using SGE/SLURM
            job_id = Chai_task.queue(
                queue='gpu.q',
                output_file=str('./job_outputs') + f'/chai_{Path(fasta_file).stem}.out',
                scheduler=scheduler,
                time='00:30:00',
                job_name=f'chai_{Path(fasta_file).stem}',
                node_name = node_name)
            print('Queueing Chai Job: ', job_id)
            chai_job_id_list.append(job_id)

    # Wait for all jobs to finish
    if scheduler == 'local':
        # If running locally, we don't need to wait for jobs.
        print('All jobs completed.')
    elif scheduler in ['SGE', 'SLURM']:
        # Wait for all jobs to finish
        print('Waiting for all jobs to finish...')
        ribbon.wait_for_jobs(chai_job_id_list, scheduler=scheduler, max_wait=60*60*24)

def fix_bonds(Chai_dir: Path, Fixed_bond_dir: Path, resnames: list, sdf_files: list):
    """
    For each structure produced by CHAI-1, fix the bonds using the provided SDF files.
    """
    for subdir in Chai_dir.iterdir():
        if subdir.is_dir():
            # Create a corresponding output directory for fixed bonds
            target_subdir = Fixed_bond_dir / subdir.stem
            make_directories(target_subdir)
            for cif_file in list_files(subdir, '.cif'):
                infile = Path(cif_file)
                for resname, sdf_file in zip(resnames, sdf_files):
                    out_file = target_subdir / Path(cif_file).name
                    patch_files(
                        cif_in=str(infile),
                        cif_out=str(out_file),
                        sdf_in=str(sdf_file),
                        resname=str(resname),
                    )
                    infile = Path(out_file) #patch multiple into same file.

def add_hydrogens(Chai_dir: Path, hydrogens_dir: Path):
    """
    Add hydrogens to each structure produced by CHAI-1.
    """
    for subdir in Chai_dir.iterdir():
        if subdir.is_dir():
            # Create a corresponding output directory for hydrogens
            target_subdir = hydrogens_dir / subdir.stem
            make_directories(target_subdir)
            for cif_file in list_files(subdir, '.cif'):
                out_file = target_subdir / Path(cif_file).name
                ribbon.AddHydrogens(
                    input_file=cif_file,
                    output_file=out_file,
                ).run()

def calculate_distance(hydrogens_dir: Path, distance_dir: Path):
    """
    Calculate the average distance (with a force restraint) between ligand and 
    specific protein atoms for each design and save the average distance.
    """
    for subdir in hydrogens_dir.iterdir():
        if subdir.is_dir():
            # Create a directory to store distance calculation results for this design
            design_distance_dir = distance_dir / subdir.stem
            make_directories(design_distance_dir)
            for cif_file in list_files(subdir, '.cif'):
                out_file_OD1 = design_distance_dir / f'{Path(cif_file).stem}.distOD1'
                out_file_OD2 = design_distance_dir / f'{Path(cif_file).stem}.distOD2'
                
                # This is slightly janky, because I've ported it backwards from my dimer code. 
                # A future monomer release will be cleaner :)
                ribbon.CalculatePairwiseDistance(
                    pdb_file=cif_file,
                    atom_list_A=['B:1:N4_1'],
                    atom_list_B=['A:129:OD1'],
                    output_file=out_file_OD1,
                    average=True,
                ).run()
                ribbon.CalculatePairwiseDistance(
                    pdb_file=cif_file,
                    atom_list_A=['B:1:N4_1'],
                    atom_list_B=['A:129:OD2'],
                    output_file=out_file_OD2,
                    average=True,
                ).run()

            # After processing, compute the average distance for the design
            distance_files_OD1 = list_files(design_distance_dir, '.distOD1')
            distances = []

            # For each pair of distance files (same stem, different suffix), get the minimum:
            for dfile_OD1 in distance_files_OD1:
                dfile_OD2 = dfile_OD1.replace('.distOD1', '.distOD2')
                with open(dfile_OD1) as f:
                    OD1_distance = float(f.read())
                with open(dfile_OD2) as f:
                    OD2_distance = float(f.read())
                # Append the minimum distance to the list
                distances.append(min(OD1_distance, OD2_distance))
            
            if distances:
                average_distance = sum(distances) / len(distances)
                with open(distance_dir / f'{subdir.stem}.avg', 'w') as avg_file:
                    avg_file.write(str(average_distance))

def get_top_designs(Chai_dir: Path, distance_dir: Path, top_design_dir: Path, top_num: int, cycle_num: int):
    """
    Identify the top designs based on the lowest average distance and copy the 
    corresponding structures to a designated directory. Also output a CSV.
    """
    # Gather average distance files.
    avg_files = list_files(distance_dir, '.avg')
    distance_dict = {}
    for file in avg_files:
        with open(file) as f:
            avg = float(f.read())
        design_name = Path(file).stem
        distance_dict[design_name] = avg

    # Sort designs by lowest average distance.
    sorted_designs = sorted(distance_dict.items(), key=lambda x: x[1])
    
    # Write results to a CSV file.
    csv_path = top_design_dir / 'top_designs.csv'
    with open(csv_path, 'w') as csv_file:
        csv_file.write('Design, Top_Name, Distance\n')
        for i, (design, dist) in enumerate(sorted_designs):
            new_name = f'cycle_{cycle_num}_top_{i}'
            csv_file.write(f'{design}, {new_name}, {dist}\n')

    # Get the top "top_num" designs.
    top_designs = dict(sorted_designs[:top_num])
    
    # For each top design, copy its predicted structure files.
    for i, design in enumerate(top_designs):
        design_dir = Chai_dir / design
        for j, cif_file in enumerate(list_files(design_dir, '.cif')):
            new_name = f'cycle_{cycle_num}_top_{i}'
            destination = top_design_dir / f'cycle_{cycle_num}_top_{i}_model_{j}.cif'
            shutil.copy(cif_file, destination)

            # Hotfix for Prody bug: reformat lines if an ATOM/HETATM line is split.
            with open(destination, 'r') as f:
                lines = f.readlines()
            with open(destination, 'w') as f:
                for line in lines:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        # If the line has too few fields (<19), join the next line.
                        if len(line.split()) < 19:
                            f.write(line.strip() + ' ')
                        else:
                            f.write(line)
                    else:
                        f.write(line)


def sculpt(input_dir: Path, run_dir: Path, num_cycles: int, input_smiles: str,
           residues_to_change: str, lmpnn_design_num: int, top_structures_num: int, scheduler: str = 'local', queue: str = '*', node_name: str = '*',
           distances: list = [], torsions: list = [], angles: list = []):
    """
    Main function to run the sculpting pipeline.

    Parameters:
    - input_dir: Directory containing the input files.
    - run_dir: Directory to store the results.
    - num_cycles: Number of cycles to run.
    - input_smiles: SMILES string of the ligand.
    - residues_to_change: Residues to be redesigned, using chain ID and residue number. Format is "A8 A39 B42 C13 C91 ...".
    - lmpnn_design_num: Number of sequences to generate per structure.
    - top_structures_num: Number of top designs to keep per cycle. Since Chai-1 generates 5 structures per design, the total structures kept per cycle is 5*top_structures_num.
    - scheduler: EasyMD and Chai-1 jobs may be submitted to a job scheduler. 'local' will run jobs locally, 'SGE' will submit jobs to the SGE scheduler, 'SLURM' will submit jobs to the SLURM scheduler.

    """
    NUM_LIGANDS = 1

    # ----------------------------
    # SETUP INITIAL DIRECTORIES
    # ----------------------------
    cycle_start_dir = run_dir / 'cycle_start'
    initial_structures_dir = cycle_start_dir / '5_Top_Designs'
    make_directories(run_dir, initial_structures_dir)
    copy_input_files(input_dir, initial_structures_dir)
    
    # Set the starting point for the first cycle.
    previous_dir = cycle_start_dir

    # ----------------------------
    # MAIN PIPELINE LOOP
    # ----------------------------
    for i in range(num_cycles):
        print('-----------------------------------')
        print('---------- Cycle:', i, '-------------')
        print('-----------------------------------')

        # Define cycle-specific directories.
        current_dir = run_dir / f'cycle_{i}'
        em_dir = current_dir / '0_EM'
        LMPNN_dir = current_dir / '1_LMPNN'
        Chai_dir = current_dir / '2_Chai'
        Fixed_bonds_dir = current_dir / '2_Chai_fixed_bonds'
        hydrogens_dir = current_dir / '3_Hydrogens'
        distance_dir = current_dir / '4_Distance'
        top_design_dir = current_dir / '5_Top_Designs'
        make_directories(current_dir, em_dir, LMPNN_dir, Chai_dir, hydrogens_dir, distance_dir, top_design_dir)

        # Step 1: Energy Minimization.
        energy_minimization(previous_dir, em_dir, '6NT_ideal.sdf', scheduler=scheduler, queue=queue, node_name=node_name, custom_bonds=distances, custom_torsions=torsions, custom_angles=angles)

        # Step 2: Run LigandMPNN for sequence design.
        run_ligand_mpnn(em_dir, LMPNN_dir, lmpnn_design_num, residues_to_change, temperature=0.2)

        # Step 3: Run CHAI-1 for structure prediction.
        run_chai(LMPNN_dir, Chai_dir, input_smiles, NUM_LIGANDS, scheduler=scheduler, queue=queue, node_name=node_name)

        # Step 3.5: Fix bonds in the generated structures.
        fix_bonds(Chai_dir, Fixed_bonds_dir, ['LIG2'], ['6NT_ideal.sdf'])

        # Step 4: Add hydrogens to predicted structures.
        add_hydrogens(Fixed_bonds_dir, hydrogens_dir)

        # Step 5: Calculate distances.
        calculate_distance(hydrogens_dir, distance_dir)

        # Step 6: Get top designs based on the distance criteria.
        get_top_designs(Chai_dir, distance_dir, top_design_dir, top_structures_num, i)

        # Prepare for the next cycle.
        previous_dir = current_dir