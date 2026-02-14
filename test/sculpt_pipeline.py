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
                     /         LIG2          /                 \                 
                 ---C                      /                   \\               
           ------    --                   /                      \\             
       O1---            --                /                         C            
       |                 --       ------C                         /             
       |                   -C-----       \\                      /              
       |                   /               \                    /               
       |                 //                 \                   /               
       |                /                    \\                /                
       N1---            /                       \        ------C                 
           ------    //                         O-------       \\               
                 ---C1                                           \\             
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
                                      ASP-129                                                  

"""
from sculpt import Sculpt, SculptOptimizer, SculptResequencer, SculptFolder
from sculpt.geometry import Atom, CustomBond, CustomAngle, CustomTorsion
from sculpt.tasks.score import SculptGeometricScoringFunction, ScoreBond, ScoreAngle
#from sculpt.sculpt import sculpt
#from sculpt.geometry import Atom, CustomBond, CustomAngle, CustomTorsion
#from sculpt.optimize import SculptOptimizer
#from sculpt.score import SculptGeometricScoringFunction, ScoreBond
from pathlib import Path


# Here are the atoms we're working with:
### First protein chain and ligand (A and X)
A_OD1 = Atom(chain='A', residue=37, name='OD1')
A_OD2 = Atom(chain='A', residue=37, name='OD2')
A_CG =  Atom(chain='A', residue=37, name='CG')
C_H =   Atom(chain='B', residue=1, name='H1')  # Ligand A
C_C3 =  Atom(chain='B', residue=1, name='C1')  # Ligand A

# And the forces we're applying to optimize:
### Get the proton close to Asp-oxygen:
A_X_bond = CustomBond(
    atom_1=A_OD1,
    atom_2=C_H,
    force_constant=20,
    target_distance=1.5
)
# Make sure the oxygen at 120 degrees, so the "bunny ears" orbital overlaps the proton:
A_X_angle = CustomAngle(
    atom_1=A_CG,
    atom_2=A_OD1,
    atom_3=C_H,
    force_constant=50,
    target_angle=120
)
# Make sure the proton's vibration vector is pointed towards the oxygen:
A_X_angle2 = CustomAngle(
    atom_1=A_OD1,
    atom_2=C_H,
    atom_3=C_C3,
    force_constant=50,
    target_angle=180
)
# And make sure the proton is in the right plane with respect to the Asp-oxygen and the CG:
A_X_torsion = CustomTorsion(
    atom_1=A_OD2,
    atom_2=A_CG,
    atom_3=A_OD1,
    atom_4=C_H,
    force_constant=20,
    periodicity=2,
    target_angle=0
)


# Create the optimizer object with the custom restraints.
optimizer = SculptOptimizer(custom_bonds=[A_X_bond],
                            custom_torsions=[A_X_torsion],
                            custom_angles=[A_X_angle, A_X_angle2])

######################################################################################
######################### DEFINING THE SCORING FUNCTION ##############################
######################################################################################

# The ScoreBonds objects let us calculate how far from the desired geometry the designed structures are.
Score_A_X_bond = ScoreBond(
    bond=A_X_bond
)

# Then, we create a scoring function where we pass in the structure.
# This function gets the average distance from the target distance for each bond.
scoring_function = SculptGeometricScoringFunction(
    function = lambda structure: (Score_A_X_bond(structure))**2 
)

########################################################################################
################################ INPUTS AND OUTPUTS ####################################
########################################################################################
input_structures_dir = Path('./data/from_run_14/')
run_dir = Path('./run_15/')

# Ligand information
ligand_reference = input_structures_dir / 'KEMP1_h.sdf'
input_smiles = 'CC1=CC(=O)OC2=C1C=CC3=C2[NH]N=N3' #KEMP1_TSA

# Which residues do we allow to change?
fixed_residues = "A37"

# Pipeline parameters
num_cycles = 15
lmpnn_design_num = 10      # Number of sequences to generate per structure (x5 Chai per sequence)
top_designs_num = 3       # Number of top designs to keep per cycle
#lmpnn_design_num = 1      # Number of sequences to generate per structure (x5 Chai per sequence)
#top_designs_num = 1       # Number of top designs to keep per cycle

structure_input_dir = Path(input_structures_dir)
current_input_designs = []

from sculpt.tasks.fold import SculptHydrogenAdder
from sculpt.design import Design

# new_design = Design(name='test')
# new_design.load_structure(Path('run_1') / 'cycle_0' / '1_Optimize' / '5RGA_KEMP1_input_opt.cif')
# new_design.structure.standardize(standard_molecules=[], renumber=True)
# new_design.structure_file('test.cif')

Sculpt( optimizer=optimizer,
       resequencer=SculptResequencer(model='LigandMPNN', num_sequences=lmpnn_design_num, fixed_residues=fixed_residues),
       folder=SculptFolder(model='Chai-1', add_hydrogens=True),
       scoring_function=scoring_function,
       structure_input_dir=input_structures_dir,
       num_cycles=num_cycles,
       run_dir=run_dir,
       sdf_files=[ligand_reference],
       top_designs_num=top_designs_num
       )