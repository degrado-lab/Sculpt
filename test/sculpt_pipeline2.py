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
from sculpt import Sculpt, SculptOptimizer, SculptResequencer, SculptFolder, DummyFolder
from sculpt.geometry import Atom, CustomBond, CustomAngle, CustomTorsion
from sculpt.tasks.score import SculptGeometricScoringFunction, ScoreBond, ScoreAngle, ScoreDihedral
#from sculpt.sculpt import sculpt
#from sculpt.geometry import Atom, CustomBond, CustomAngle, CustomTorsion
#from sculpt.optimize import SculptOptimizer
#from sculpt.score import SculptGeometricScoringFunction, ScoreBond
from pathlib import Path

import numpy as np


# Here are the atoms we're working with:
### First protein chain and ligand (A and X)
A_OD1 = Atom(chain='A', residue=99, name='OD1')
A_OD2 = Atom(chain='A', residue=99, name='OD2')
A_CG =  Atom(chain='A', residue=99, name='CG')
B_H =   Atom(chain='B', residue=1, name='H7')  # Ligand A
B_N3 =  Atom(chain='B', residue=1, name='N3')  # Ligand A

# And the forces we're applying to optimize:
### Get the proton close to Asp-oxygen:
A_X_bond = CustomBond(
    atom_1=A_OD1,
    atom_2=B_H,
    force_constant=20,
    target_distance=1.16
)
# Make sure the oxygen at 120 degrees, so the "bunny ears" orbital overlaps the proton:
A_X_angle = CustomAngle(
    atom_1=A_CG,
    atom_2=A_OD1,
    atom_3=B_H,
    force_constant=150, #Upped from 50
    target_angle=120
)
# Make sure the proton's vibration vector is pointed towards the oxygen:
A_X_angle2 = CustomAngle(
    atom_1=A_OD1,
    atom_2=B_H,
    atom_3=B_N3,
    force_constant=150, # Upped from 50
    target_angle=180
)
# And make sure the proton is in the right plane with respect to the Asp-oxygen and the CG:
### BUG! In EasyMD, a torsion angle of 180 goes to 0 (and 0 goes to 180).
### So to get to 0, we'll define a 180 target angle here, and then in scoring, we'll score against 0.
A_X_torsion = CustomTorsion(
    atom_1=A_OD2,
    atom_2=A_CG,
    atom_3=A_OD1,
    atom_4=B_H,
    force_constant=180, #Upped from 80
    periodicity=1,
    target_angle=180
)
A_X_torsion_for_scoring = CustomTorsion(
    atom_1=A_OD2,
    atom_2=A_CG,
    atom_3=A_OD1,
    atom_4=B_H,
    force_constant=180, #Upped from 80
    periodicity=1,
    target_angle=0
)


# Create the optimizer object with the custom restraints.
optimizer = SculptOptimizer(custom_bonds=[A_X_bond],
                            custom_torsions=[A_X_torsion],
                            custom_angles=[A_X_angle, A_X_angle2])
                            #ull_sim=True)  # Since we have difficult constraints, let's shake out the protein.

######################################################################################
######################### DEFINING THE SCORING FUNCTION ##############################
######################################################################################

# The ScoreBonds objects let us calculate how far from the desired geometry the designed structures are.
Score_A_X_bond = ScoreBond(
    bond=A_X_bond
)

Score_A_X_angle = ScoreAngle(
    angle=A_X_angle,
    absolute=True
)
Score_A_X_angle2 = ScoreAngle(
    angle=A_X_angle2,
    absolute=True
)

Score_A_X_dihedral = ScoreDihedral(
    torsion=A_X_torsion_for_scoring,  # applying our bug-patch here
    absolute=True
)

# Then, we create a scoring function where we pass in the structure.
# This function gets the average distance from the target distance for each bond.

def cos_adjust(x):
    return (1/2) * (1 + np.cos( (2*np.pi * (x))/360 )) # I used to subtract 180 from x, but this keeps the best score = 1 instead of 0.

scoring_function = SculptGeometricScoringFunction(
    function = lambda design: max(20-Score_A_X_bond(design), 0) + cos_adjust(Score_A_X_angle(design)) + cos_adjust(Score_A_X_angle2(design)) + cos_adjust(Score_A_X_dihedral(design))
) # Squared the bond score to make it more important (outweighs the angles better)



########################################################################################
################################ INPUTS AND OUTPUTS ####################################
########################################################################################

input_structures_dir = Path('./data/1OHP_monomer/')
run_dir = Path('./1OHP_no_opt_30pop/')

# Ligand information
ligand_reference = input_structures_dir / 'KEMP1_TSA_h.sdf'

# Which residues do we allow to change?
fixed_residues = "A99"

# Pipeline parameters
num_cycles = 15
lmpnn_design_num = 10      # Number of sequences to generate per structure (x5 Chai per sequence)
top_designs_num = 3       # Number of top designs to keep per cycle

structure_input_dir = Path(input_structures_dir)
current_input_designs = []

from sculpt.tasks.fold import SculptHydrogenAdder
from sculpt.design import Design

# from sculpt.design import Design
# new_design = Design(name='test')
# new_design.load_structure(Path('data') / 'miscreant_cordon_HHH18' / 'Miscreant_Cordon_b86ca4c9_fold_0.cif')
# new_design.structure.standardize(standard_molecules=[], renumber=True)
# new_design.structure_file('test.cif')

# # Test scoring function:
# print('Score components:')
# a =  max(10-Score_A_X_bond(new_design), 0)
# b = cos_adjust(Score_A_X_angle(new_design))
# c = cos_adjust(Score_A_X_angle2(new_design))
# d = cos_adjust(Score_A_X_dihedral(new_design))
# print(a, b, c, d)
# print('Total score:', a + b + c + d)
# #print(scoring_function.score(new_design))

Sculpt( optimizer=optimizer,
       resequencer=SculptResequencer(model='LASErMPNN', num_sequences=lmpnn_design_num, fixed_residues=fixed_residues),
       folder=SculptFolder(model='Chai-1', add_hydrogens=True),
       #folder=DummyFolder(),
       scoring_function=scoring_function,
       # Random number 0 to 10:
       #scoring_function=SculptGeometricScoringFunction(function = lambda design: np.random.uniform(0, 10)),
       structure_input_dir=input_structures_dir,
       num_cycles=num_cycles,
       run_dir=run_dir,
       sdf_files=[ligand_reference],
       top_designs_num=top_designs_num,
       cycle_retry=5,
       pop_size=30,
       )

