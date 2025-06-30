def sculpt(input_structures_dir, run_dir, ligand_reference, input_smiles, num_ligands,
           residues_to_change, num_cycles, lmpnn_design_num, top_designs_num, optimizer, scoring_function):
    """Main function to run the Sculpt pipeline for ligand design and optimization.
    
    This function orchestrates the complete Sculpt pipeline, including energy minimization,
    sequence design with LigandMPNN, structure prediction with CHAI-1, hydrogen addition,
    scoring, and iterative refinement over multiple cycles.
    
    Args:
        input_structures_dir: Directory containing input structures.
        run_dir: Directory to store the results.
        ligand_reference: Reference ligand SDF file.
        input_smiles: SMILES string of the input ligand.
        num_ligands: Number of ligands to design.
        residues_to_change: Residues to allow for redesign.
        num_cycles: Number of cycles to run the design pipeline.
        lmpnn_design_num: Number of sequences to generate per structure.
        top_designs_num: Number of top designs to keep per cycle.
        optimizer: Optimizer object for structure optimization.
        scoring_function: Scoring function for evaluating designs.
    """
    # -----------------------------
    # CONSTANTS & INPUT PARAMETERS
    # -----------------------------
    
    # Ligand information
    ligand_reference = 'KEMP1_TSA_h.sdf'
    INPUT_SMILES = input_smiles
    NUM_LIGANDS = num_ligands

    # Protein redesign information: keep your interface residues constant.
    RESIDUES_TO_CHANGE = residues_to_change

    # Pipeline parameters
    NUM_CYCLES = num_cycles
    LMPNN_DESIGN_NUM = lmpnn_design_num      # Number of sequences to generate per structure
    TOP_STRUCTURES_NUM = top_designs_num     # Number of top designs to keep per cycle

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


######################################################################################
####################### DEFINING THE GEOMETRY TO OPTIMIZE ############################
######################################################################################

"""
                                                                             
                                                                C               
                                                               /                
                         C------------C                       /                 
                        /              \\                    /                  
                        /                \                   /                  
                       /                  \                 /                   
                      /                    \\       -------C                    
                     /                       C------        \\                  
                     /         LIG          /                 \                 
                 ---C                      /                   \\               
           ------    --                   /                      \\             
       N---            --                /                         C            
       |                 --       ------C                         /             
       |                   -C-----       \\                      /              
       |                   /               \                    /               
       |                 //                 \                   /               
       |                /                    \\                /                
       N---            /                       \        ------C                 
           ------    //                         O-------       \\               
                 ---N3                                           \\             
                    |                                              \            
                    |                                               \\          
                    |                 OD2                             O         
                    H7                /                                       
                                     /                                         
                                    /
                    OD1-----------CG                                       
                                    \                                
                                     \                                         
                                      \                                       
                                      ASP-99                                                  

"""
# Here are the atoms we're working with:
### First protein chain and ligand (A and X)
A_OD1 = Atom(chain='A', residue=99, name='OD1')
A_OD2 = Atom(chain='A', residue=99, name='OD2')
A_CG =  Atom(chain='A', residue=99, name='CG')
X_H =   Atom(chain='X', residue=1, name='H7')  # Ligand A
X_N3 =  Atom(chain='X', residue=1, name='N3')  # Ligand A

### Second protein chain and ligand (B and Y)
B_OD1 = Atom(chain='B', residue=99, name='OD1')
B_OD2 = Atom(chain='B', residue=99, name='OD2')
B_CG =  Atom(chain='B', residue=99, name='CG')
Y_H =   Atom(chain='Y', residue=1, name='H7')  # Ligand B
Y_N3 =  Atom(chain='Y', residue=1, name='N3')  # Ligand B

# And the forces we're applying to optimize:
### Get the proton close to Asp-oxygen:
A_X_bond = CustomBond(
    atom_1=A_OD1,
    atom_2=X_H,
    force_constant=20,
    target_distance=1.5
)
# Make sure the oxygen at 120 degrees, so the "bunny ears" orbital overlaps the proton:
A_X_angle = CustomAngle(
    atom_1=A_CG,
    atom_2=A_OD1,
    atom_3=X_H,
    force_constant=50,
    target_angle=120
)
# Make sure the proton's vibration vector is pointed towards the oxygen:
A_X_angle2 = CustomAngle(
    atom_1=A_OD1,
    atom_2=X_H,
    atom_3=X_N3,
    force_constant=50,
    target_angle=180
)
# And make sure the proton is in the right plane with respect to the Asp-oxygen and the CG:
A_X_torsion = CustomTorsion(
    atom_1=A_OD2,
    atom_2=A_CG,
    atom_3=A_OD1,
    atom_4=X_H,
    force_constant=20,
    periodicity=2,
    target_angle=0
)

### Repeat the exact same geometry for the second ligand and second protein chain.
B_Y_bond = CustomBond(
    atom_1=B_OD1,
    atom_2=Y_H,
    force_constant=20,
    target_distance=1.5
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
B_Y_torsion = CustomTorsion(
    atom_1=B_OD2,
    atom_2=B_CG,
    atom_3=B_OD1,
    atom_4=Y_H,
    force_constant=20,
    periodicity=2,
    target_angle=0
)


# Create the optimizer object with the custom restraints.
optimizer = SculptOptimizer(custom_bonds=[A_X_bond, B_Y_bond],
                            custom_torsions=[A_X_torsion, B_Y_torsion],
                            custom_angles=[A_X_angle, A_X_angle2, B_Y_angle, B_Y_angle2])

######################################################################################
######################### DEFINING THE SCORING FUNCTION ##############################
######################################################################################

# The ScoreBonds objects let us calculate how far from the desired geometry the designed structures are.
Score_A_X_bond = ScoreBond(
    bond=A_X_bond
)
Score_B_Y_bond = ScoreBond(
    bond=B_Y_bond
)

# Then, we create a scoring function where we pass in the structure.
# This function gets the average distance from the target distance for each bond.
scoring_function = SculptGeometricScoringFunction(
    function = lambda structure: (Score_A_X_bond(structure) + Score_B_Y_bond(structure)) / 2
)

########################################################################################
################################ INPUTS AND OUTPUTS ####################################
########################################################################################
input_structures_dir = Path('./inputs/')
run_dir = Path('./run_1/')

# Ligand information
ligand_reference = 'KEMP1_TSA_h.sdf'
input_smiles = 'CC1=CC(=O)OC2=C1C=CC3=C2[NH]N=N3' #KEMP1_TSA
num_ligands = 2

# Which residues do we allow to change?
residues_to_change = (
    "A14 A18 A38 A54 A55 A58 A59 A63 A82 A84 A86 A95 A112 A114 A115 A116" +
    "B14 B18 B38 B54 B55 B58 B59 B63 B82 B84 B86 B95 B112 B114 B115 B116"
)

# Pipeline parameters
num_cycles = 15
lmpnn_design_num = 1      # Number of sequences to generate per structure
top_designs_num = 1     # Number of top designs to keep per cycle

sculpt()

if __name__ == '__main__':
    main()