import json
import shutil
from pathlib import Path
import ribbon
from ribbon.utils import make_directories, list_files
import shim
import shim.convert
from sculpt_helpers import SculptOptimizer, Atom, CustomBond, CustomAngle, CustomTorsion, ScoreBond, SculptGeometricScoringFunction


def copy_input_files(input_dir: Path, target_dir: Path):
    """Copy all files from the input directory to the target directory.
    
    Args:
        input_dir: Source directory containing files to copy.
        target_dir: Destination directory for the copied files.
    """
    for file in input_dir.iterdir():
        shutil.copy(file, target_dir)


def energy_minimization(previous_dir: Path, em_dir: Path, sdf_file: str, optimizer: SculptOptimizer):
    """Energy minimization for all structures from the previous cycle.
    
    For every structure in the previous cycle, run energy minimization using
    custom restraints defined in the optimizer.
    
    Args:
        previous_dir: Directory containing structures from the previous cycle.
        em_dir: Directory to store energy-minimized structures.
        sdf_file: Path to the ligand SDF file for parameterization.
        optimizer: SculptOptimizer object containing custom restraints.
    """
    # Look for both .cif and .pdb files in the previous design directory.
    prev_design_dir = previous_dir / '5_Top_Designs'
    structure_files = list_files(prev_design_dir, '.cif') + list_files(prev_design_dir, '.pdb')
    print('Number of structures:', len(structure_files))
    
    for file in structure_files:
        output_prefix = em_dir / Path(file).stem
        ligA = "X" # based on our A->X, B->Y renaming scheme
        ligB = "Y"
        
        # Run energy minimization with custom restraints
        ribbon.EasyMD(
            input_file=file,
            output_prefix=output_prefix,
            ligand_files=[sdf_file],
            duration=1,
            custom_bonds=optimizer.get_custom_bonds(),
            custom_torsions=optimizer.get_custom_torsions(),
            custom_angles=optimizer.get_custom_angles(),
            minimize_only=True,
        ).run()


def run_ligand_mpnn(em_dir: Path, LMPNN_dir: Path, design_num: int, residues_to_change: str, temperature: float = 0.1):
    """Run LigandMPNN to generate sequences for each energy-minimized structure.
    
    Also fixes the filenames in the 'seqs_split' folder for downstream processing.
    
    Args:
        em_dir: Directory containing energy-minimized structures.
        LMPNN_dir: Directory to store LigandMPNN outputs.
        design_num: Number of sequence designs to generate per structure.
        residues_to_change: String specifying which residues to allow for redesign.
        temperature: Sampling temperature for sequence generation. Defaults to 0.1.
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


def run_chai(LMPNN_dir: Path, Chai_dir: Path, smiles: str, num_ligands: int):
    """Run CHAI-1 structure prediction for each designed sequence.
    
    For each designed sequence file, run CHAI-1 to predict protein-ligand
    complex structures.
    
    Args:
        LMPNN_dir: Directory containing LigandMPNN sequence outputs.
        Chai_dir: Directory to store CHAI-1 structure predictions.
        smiles: SMILES string of the ligand.
        num_ligands: Number of ligand copies in each complex.
    """
    seqs_split_dir = LMPNN_dir / 'seqs_split'
    for fasta_file in list_files(seqs_split_dir, '.fasta'):
        fasta_file = Path(fasta_file)
        ribbon.Chai1(
            fasta_file=fasta_file,
            smiles_string=smiles,
            num_ligands=num_ligands,
            output_dir=Chai_dir / fasta_file.stem,
        ).run()


def add_hydrogens(Chai_dir: Path, hydrogens_dir: Path, db_file: Path):
    """Add hydrogens to each structure produced by CHAI-1.
    
    Uses the Reduce tool to add hydrogen atoms to all PDB structures
    in the CHAI-1 output directories.
    
    Args:
        Chai_dir: Directory containing CHAI-1 structure predictions.
        hydrogens_dir: Directory to store structures with added hydrogens.
        db_file: Database file for Reduce containing ligand information.
    """
    for subdir in Chai_dir.iterdir():
        if subdir.is_dir():
            # Create a corresponding output directory for hydrogens
            target_subdir = hydrogens_dir / subdir.stem
            make_directories(target_subdir)
            for pdb_file in list_files(subdir, '.pdb'):
                out_file = target_subdir / Path(pdb_file).name
                ribbon.Reduce(
                    pdb_input_file=pdb_file,
                    pdb_output_file=out_file,
                    database=db_file,
                ).run()


def score(hydrogens_dir: Path, distance_dir: Path, scoring_function: callable):
    """Score structures using the provided scoring function.
    
    Applies the scoring function to each structure and saves the results
    to score files for downstream analysis.
    
    Args:
        hydrogens_dir: Directory containing structures with hydrogens added.
        distance_dir: Directory to store scoring results.
        scoring_function: Callable function that takes a structure file and returns a score.
    """
    for subdir in hydrogens_dir.iterdir():
        if subdir.is_dir():
            # Create a directory to store distance calculation results for this design
            #design_distance_dir = distance_dir / subdir.stem
            #make_directories(design_distance_dir)
            score_list = []
            for cif_file in list_files(subdir, '.pdb'):
                # out_file = design_distance_dir / f'{Path(cif_file).stem}.dist'
                # ribbon.CalculatePairwiseDistance(
                #     pdb_file=cif_file,
                #     atom_list_A=['X:1:H7', 'Y:1:H7'],
                #     atom_list_B=['A:99:OD1', 'A:99:OD2', 'B:99:OD1', 'B:99:OD2'],
                #     output_file=out_file,
                #     average=True,
                # ).run()
                current_score = scoring_function(cif_file)
                score_list.append(current_score)

            # After processing, compute the average distance for the design
            #distance_files = list_files(design_distance_dir, '.dist')
            #distances = []
            #for dfile in distance_files:
            #    with open(dfile) as f:
            #        distances.append(float(f.read()))
            if score_list:
                average_score = sum(score_list) / len(score_list)
                with open(distance_dir / f'{subdir.stem}.avg', 'w') as avg_file:
                    avg_file.write(str(average_score))

def calculate_distance(hydrogens_dir: Path, distance_dir: Path):
    """Calculate average distances between ligand and protein atoms.
    
    Calculate the average distance (with a force restraint) between ligand and 
    specific protein atoms for each design and save the average distance.
    
    Args:
        hydrogens_dir: Directory containing structures with hydrogens added.
        distance_dir: Directory to store distance calculation results.
    """
    for subdir in hydrogens_dir.iterdir():
        if subdir.is_dir():
            # Create a directory to store distance calculation results for this design
            design_distance_dir = distance_dir / subdir.stem
            make_directories(design_distance_dir)
            for cif_file in list_files(subdir, '.pdb'):
                out_file = design_distance_dir / f'{Path(cif_file).stem}.dist'
                ribbon.CalculatePairwiseDistance(
                    pdb_file=cif_file,
                    atom_list_A=['X:1:H7', 'Y:1:H7'],
                    atom_list_B=['A:99:OD1', 'A:99:OD2', 'B:99:OD1', 'B:99:OD2'],
                    output_file=out_file,
                    average=True,
                ).run()

            # After processing, compute the average distance for the design
            distance_files = list_files(design_distance_dir, '.dist')
            distances = []
            for dfile in distance_files:
                with open(dfile) as f:
                    distances.append(float(f.read()))
            if distances:
                average_distance = sum(distances) / len(distances)
                with open(distance_dir / f'{subdir.stem}.avg', 'w') as avg_file:
                    avg_file.write(str(average_distance))


def get_top_designs(Chai_dir: Path, distance_dir: Path, top_design_dir: Path, top_num: int, cycle_num: int):
    """Identify and copy the top designs based on scoring metrics.
    
    Identify the top designs based on the lowest average distance and copy the 
    corresponding structures to a designated directory. Also output a CSV summary.
    
    Args:
        Chai_dir: Directory containing CHAI-1 structure predictions.
        distance_dir: Directory containing scoring/distance results.
        top_design_dir: Directory to store the top-ranked designs.
        top_num: Number of top designs to select.
        cycle_num: Current cycle number for logging purposes.
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
        for j, pdb_file in enumerate(list_files(design_dir, '.pdb')):
            new_name = f'cycle_{cycle_num}_top_{i}'
            destination = top_design_dir / f'cycle_{cycle_num}_top_{i}_model_{j}.pdb'
            shutil.copy(pdb_file, destination)

            # Hotfix for Prody bug: reformat lines if an ATOM/HETATM line is split.
            # with open(destination, 'r') as f:
            #     lines = f.readlines()
            # with open(destination, 'w') as f:
            #     for line in lines:
            #         if line.startswith('ATOM') or line.startswith('HETATM'):
            #             # If the line has too few fields (<19), join the next line.
            #             if len(line.split()) < 19:
            #                 f.write(line.strip() + ' ')
            #             else:
            #                 f.write(line)
            #         else:
            #             f.write(line)

def standardize_structure_dir(structure_dir: Path, resname: str, ligand_reference: str, ligand_chain_mapping = {'A':'X', 'B':'Y'}, use_hydrogens=True, rename_from=None):
    """Rename ligand atoms in structures to match the reference ligand.
    
    Processes all structure files in a directory to standardize ligand atom names
    and chain identifiers according to the reference ligand.
    
    Args:
        structure_dir: Directory containing structure files to process.
        resname: Residue name for the ligand.
        ligand_reference: Path to the reference ligand file.
        ligand_chain_mapping: Mapping of protein chains to ligand chains. Defaults to {'A':'X', 'B':'Y'}.
        use_hydrogens: Whether to include hydrogen atoms. Defaults to True.
        rename_from: List of alternative residue names to rename from. Defaults to None.
    """
    for file in list_files(structure_dir, '.cif') + list_files(structure_dir, '.pdb'):
        print('Renaming ligand atoms in:', file)
        standardize_structure(file, file, resname, ligand_reference, ligand_chain_mapping, use_hydrogens=use_hydrogens, rename_from=rename_from)

def standardize_structure(structure_file: Path, output_file: Path, resname: str, ligand_reference: str, ligand_chain_mapping = {'A':'X', 'B':'Y'}, use_hydrogens=True, rename_from= None):
    """Standardize a single structure file.
    
    Standardize the structure by renaming ligand atoms and chains to match the reference ligand.
    
    Args:
        structure_file: Input structure file to standardize.
        output_file: Output file for the standardized structure.
        resname: Residue name for the ligand.
        ligand_reference: Path to the reference ligand file.
        ligand_chain_mapping: Mapping of protein chains to ligand chains. Defaults to {'A':'X', 'B':'Y'}.
        use_hydrogens: Whether to include hydrogen atoms. Defaults to True.
        rename_from: List of alternative residue names to rename from. Defaults to None.
    """
    ## Rename all lig residues to LIG:
    if rename_from is not None:
        extension = Path(structure_file).suffix.lower()
        for old_resname in rename_from:
            in_structure = shim.ShimStructure(cif_file = structure_file) if extension == '.cif' else shim.ShimStructure(pdb_file = structure_file)
            in_structure.rename_residue(old_resname, resname)
            in_structure.to_cif(structure_file)

    # First, rename ligand atoms to match the reference ligand.
    shim.fix.fix_atom_names_in_residue(
        infile=structure_file,
        outfile=output_file,
        resname=resname,
        sdf_file=ligand_reference,
        use_hydrogens=use_hydrogens,
    )
    
    shim.fix.rename_lig_chain_by_proximity(
        infile = output_file,
        outfile = output_file,
        resname = resname,
        chain_mapping=ligand_chain_mapping,
    )

def write_reduce_db(input_sdf: Path, output_file: Path):
    """Create a Reduce database file from a ligand SDF file.
    
    Converts a ligand SDF file into a Reduce-compatible database format
    for hydrogen addition calculations.
    
    Args:
        input_sdf: Path to the input SDF file containing the ligand.
        output_file: Path for the output Reduce database file.
    """
    # Create a StandardMolecule:

    sm = shim.StandardMolecule(sdf_file=input_sdf)

    # Write the molecule to a file in Reduce database format.
    shim.convert.write_reduce_db(standard_mol=sm,
                      residue_name='LIG',
                      chemical_name='KEMP1_TSA',
                      output_path=output_file)
    
def main():
    """Main function to run the Sculpt pipeline with predefined parameters.
    
    Sets up and executes the complete Sculpt pipeline for ligand design
    and optimization using hardcoded parameters for a specific example.
    """
    # -----------------------------
    # CONSTANTS & INPUT PARAMETERS
    # -----------------------------
    input_dir = Path('./inputs_ligandoutTSA/')
    run_dir = Path('./test8/')
    
    # Ligand information
    ligand_reference = 'KEMP1_TSA_h.sdf'
    INPUT_SMILES = 'CC1=CC(=O)OC2=C1C=CC3=C2[NH]N=N3' #KEMP1_TSA
    NUM_LIGANDS = 2

    # Protein redesign information: keep your interface residues constant.
    RESIDUES_TO_CHANGE = (
        "A14 A18 A38 A54 A55 A58 A59 A63 A82 A84 A86 A95 A112 A114 A115 A116" +
        "B14 B18 B38 B54 B55 B58 B59 B63 B82 B84 B86 B95 B112 B114 B115 B116"
    )

    # Pipeline parameters
    NUM_CYCLES = 5
    LMPNN_DESIGN_NUM = 1      # Number of sequences to generate per structure
    TOP_STRUCTURES_NUM = 1     # Number of top designs to keep per cycle

    # ----------------------------
    # SETUP INITIAL DIRECTORIES
    # ----------------------------
    cycle_start_dir = run_dir / 'cycle_start'
    initial_structures_dir = cycle_start_dir / '5_Top_Designs'
    make_directories(run_dir, initial_structures_dir)
    copy_input_files(input_dir, initial_structures_dir)

    ## New: Shim the intial structure to make atom names and residue names standard.
    standardize_structure_dir(initial_structures_dir, 'LIG', ligand_reference, use_hydrogens=True)

    # Set the starting point for the first cycle.
    previous_dir = cycle_start_dir

    # For later, create the DB object for REDUCE:
    reduce_db_file = run_dir / 'reduce_db_LIG.txt'
    write_reduce_db(
        input_sdf=Path(ligand_reference),
        output_file=reduce_db_file
    )

    # ----------------------------
    # DEFINE OPTIMIZER
    # ----------------------------
    
    # Here are the atoms we're working with:
    A_OD1 = Atom(chain='A', residue=99, name='OD1')
    A_OD2 = Atom(chain='A', residue=99, name='OD2')
    A_CG =  Atom(chain='A', residue=99, name='CG')
    X_H =   Atom(chain='X', residue=1, name='H7')  # Ligand A
    X_N3 =  Atom(chain='X', residue=1, name='N3')  # Ligand A

    B_OD1 = Atom(chain='B', residue=99, name='OD1')
    B_OD2 = Atom(chain='B', residue=99, name='OD2')
    B_CG =  Atom(chain='B', residue=99, name='CG')
    Y_H =   Atom(chain='Y', residue=1, name='H7')  # Ligand B
    Y_N3 =  Atom(chain='Y', residue=1, name='N3')  # Ligand B

    # And the forces we're applying to optimize:
    A_X_bond = CustomBond(
        atom_1=A_OD1,
        atom_2=X_H,
        force_constant=20,
        target_distance=1.5
    )
    B_Y_bond = CustomBond(
        atom_1=B_OD1,
        atom_2=Y_H,
        force_constant=20,
        target_distance=1.5
    )
    A_X_torsion = CustomTorsion(
        atom_1=A_OD2,
        atom_2=A_CG,
        atom_3=A_OD1,
        atom_4=X_H,
        force_constant=20,
        periodicity=2,
        target_angle=0
    )
    B_Y_torsion = CustomTorsion(
        atom_1=B_OD2,
        atom_2=B_CG,
        atom_3=B_OD1,
        atom_4=Y_H,
        force_constant=20,
        periodicity=2,
        target_angle=0
    )
    A_X_angle = CustomAngle(
        atom_1=A_CG,
        atom_2=A_OD1,
        atom_3=X_H,
        force_constant=50,
        target_angle=120
    )
    A_X_angle2 = CustomAngle(
        atom_1=A_OD1,
        atom_2=X_H,
        atom_3=X_N3,
        force_constant=50,
        target_angle=180
    )
    B_Y_angle = CustomAngle(
        atom_1=B_CG,
        atom_2=B_OD1,
        atom_3=Y_H,
        force_constant=50,
        target_angle=120
    )
    B_Y_angle2 = CustomAngle(
        atom_1=B_OD1,
        atom_2=Y_H,
        atom_3=Y_N3,
        force_constant=50,
        target_angle=180
    )
    # Create the optimizer object with the custom restraints.
    optimizer = SculptOptimizer(custom_bonds=[A_X_bond, B_Y_bond],
                                custom_torsions=[A_X_torsion, B_Y_torsion],
                                custom_angles=[A_X_angle, A_X_angle2, B_Y_angle, B_Y_angle2])

    print(optimizer.get_custom_bonds())

    # Create the scoring function, using the same geometric features, and up-weighting the distance terms.
    A_X_scorebond = ScoreBond(
        bond=A_X_bond
    )
    B_Y_scorebond = ScoreBond(
        bond=B_Y_bond
    )
    scoring_function = SculptGeometricScoringFunction(
        function = lambda s: A_X_scorebond.calc(s) + B_Y_scorebond.calc(s),
    )

    # ----------------------------
    # MAIN PIPELINE LOOP
    # ----------------------------
    for i in range(NUM_CYCLES):
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
        score_dir = current_dir / '4_Distance'
        top_design_dir = current_dir / '5_Top_Designs'
        make_directories(current_dir, em_dir, LMPNN_dir, Chai_dir, hydrogens_dir, score_dir, top_design_dir)

        # Step 1: Energy Minimization.
        energy_minimization(previous_dir, em_dir, ligand_reference, optimizer)

        # Step 2: Run LigandMPNN for sequence design.
        run_ligand_mpnn(em_dir, LMPNN_dir, LMPNN_DESIGN_NUM, RESIDUES_TO_CHANGE, temperature=0.2)

        # Step 3: Run CHAI-1 for structure prediction.
        run_chai(LMPNN_dir, Chai_dir, INPUT_SMILES, NUM_LIGANDS)

        # Step 3.5: Fix bonds in the generated structures.
        ## NB. This will error out if both ligands are attached to the same protein chain (we're not doing fancy bipartite matching here)
        for subdir in Chai_dir.iterdir():
            print('Processing subdir:', subdir)
            if subdir.is_dir():
                standardize_structure_dir(subdir, 'LIG', ligand_reference, use_hydrogens=False, rename_from=['LIG3', 'LIG4'])
        ## fix_bonds(Chai_dir, Fixed_bonds_dir, ['LIG3', 'LIG4'], ['KEMP1_TSA_h_fromChai_kek.sdf', 'KEMP1_TSA_h_fromChai_kek.sdf'])
        
        # Step 3.55: Convert to PDB for reduce:
        for subdir in Chai_dir.iterdir():
            if subdir.is_dir():
                for cif_file in list_files(subdir, '.cif'):
                    # Convert .cif to .pdb
                    structure = shim.ShimStructure(cif_file=cif_file)
                    pdb_file = subdir / (Path(cif_file).stem + '.pdb')
                    structure.to_pdb(str(pdb_file))

        # Step 4: Add hydrogens to predicted structures.
        add_hydrogens(Chai_dir, hydrogens_dir, reduce_db_file)

        # Step 5: Calculate distances.
        #calculate_distance(hydrogens_dir, distance_dir)
        score(hydrogens_dir, score_dir, scoring_function)

        # Step 6: Get top designs based on the distance criteria.
        get_top_designs(hydrogens_dir, score_dir, top_design_dir, TOP_STRUCTURES_NUM, i)

        # Prepare for the next cycle.
        previous_dir = current_dir


if __name__ == '__main__':
    main()
