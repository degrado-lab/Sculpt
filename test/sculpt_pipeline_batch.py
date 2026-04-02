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
from sculpt.ga.mutate import SimpleAndDirectedMutation
from sculpt import Sculpt, SculptOptimizer, SculptResequencer, SculptFolder, DummyFolder
from sculpt.geometry import Atom, CustomBond, CustomAngle, CustomTorsion
from sculpt.tasks.score import SculptGeometricScoringFunction, ScoreBond, ScoreAngle, ScoreDihedral
from sculpt.ga.select import TournamentSelection
from sculpt.ga.crossover import SpatialCrossover
from sculpt.ga.mutate import DirectedMutation, SimpleMutation
from sculpt.alter import SculptResidueFlipper
from pathlib import Path

import numpy as np

######################################################################################
####################### DEFINING OUR ATOMS ###########################################
######################################################################################

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
    target_angle=0 # flipped here! Used to be 180
)
A_X_torsion_for_scoring = CustomTorsion(
    atom_1=A_OD2,
    atom_2=A_CG,
    atom_3=A_OD1,
    atom_4=B_H,
    force_constant=180, #Upped from 80
    periodicity=1,
    target_angle=180 # flipped here! Used to be 0
)
A_X_bond_for_scoring = CustomBond(
    atom_1=A_OD1,
    atom_2=B_H,
    force_constant=20,
    target_distance=0 # New bond for scoring, lower the better!
)

######################################################################################
######################### DEFINING THE SCORING FUNCTION ##############################
######################################################################################

# The ScoreBonds objects let us calculate how far from the desired geometry the designed structures are.
Score_A_X_bond = ScoreBond(
    bond=A_X_bond_for_scoring
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

def cos_adjust_min_at_180(x): # When x = 0, score = 1. When x = 180, score = 0.
    return (1/2) * (1 + np.cos( (2*np.pi * (x))/360 ))

def cos_adjust_min_at_0(x): # When x = 180, score = 1. When x = 0, score = 0.
    return (1/2) * (1 + np.cos( (2*np.pi * (x - 180))/360 ))

# I've flipped this function to return to the "cost" model, where we want to minimize sum of distances and angles. 
# Then, at the end we subtract the cost from 20 to get a fitness score (where 20 is perfect fit).
scoring_function = SculptGeometricScoringFunction(
    scorers={
        'bond_1': Score_A_X_bond,
        'angle_1': Score_A_X_angle,
        'angle_2': Score_A_X_angle2,
        'dihedral': Score_A_X_dihedral
    },
    aggregator=lambda results: max(100 - (results['bond_1']**2 + cos_adjust_min_at_0(results['angle_1']) + cos_adjust_min_at_0(results['angle_2']) + cos_adjust_min_at_0(results['dihedral']) ), 0)
) 

######################################################################################
################# DEFINING THE SELECTION, MUTATION, and CROSSOVER OPERATORS ##########
######################################################################################

# These will be used for folding and fixing structures before scoring:
folder = SculptFolder(
        model='Boltz-2',
        num_structures=1, 
        add_hydrogens=True,
        use_queue=True,
        scheduler="SLURM"
    )
#folder = DummyFolder(data_dir='../run_data/1OHP_mutants/')

fixers = [SculptResidueFlipper(target_atom=B_H, flip_atom=A_OD1)]

# Selection:
selector = TournamentSelection(tournament_size=7) # raised from 3. Should be quite strong selection

# Mutation:
fixed_residues = "A99" # Which residues do we allow to change?
input_structures_dir = Path('./data/1OHP_monomer/')
ligand_reference = input_structures_dir / 'KEMP1_TSA_h.sdf'
# mutator = DirectedMutation(
#             mutation_rate=0.25,
#             resequencer=    SculptResequencer(model='LASErMPNN', num_sequences=1, 
#                                             fixed_residues=fixed_residues),
#             optimizer=      SculptOptimizer(custom_bonds=[A_X_bond],
#                                             custom_torsions=[A_X_torsion],
#                                             custom_angles=[A_X_angle, A_X_angle2],
#                                             full_sim=True), 
#             sdf_files=[ligand_reference],
#             folder=         folder
# )
mutator = SimpleAndDirectedMutation( 

        # Simple Mutation params
        simple_mutation_rate = 0.5,
        mutations_per_sequence = 1,
        fixed_residues = fixed_residues,
        
        # Directed Mutation params
        directed_mutation_rate = 0.25,

        resequencer =   SculptResequencer(model='LASErMPNN', num_sequences=1, 
                                        fixed_residues=fixed_residues,
                                        scheduler="SLURM",
                                        use_queue=True),
                                        
        optimizer =     SculptOptimizer(custom_bonds=[A_X_bond],
                                        custom_torsions=[A_X_torsion],
                                        custom_angles=[A_X_angle, A_X_angle2],
                                        full_sim=True, 
                                        scheduler="SLURM",
                                        use_queue=True),

        sdf_files=[ligand_reference],
        folder=         folder,
)

# Crossover:
crossover = SpatialCrossover(crossover_rate=0.5)

# This is just to fill out our population:
initial_population_mutator = SimpleMutation(mutation_rate=1.0, mutations_per_sequence=30, temperature=2.0, fixed_residues=fixed_residues)


########################################################################################
################################ INPUTS AND OUTPUTS ####################################
########################################################################################

input_structures_dir = Path('./data/1OHP_monomer/')
run_dir = Path('./test4_1OHP_combinedmutants25p_spatialxover_30pop_dsquared100max_7xtournament_fix_boltz2_laser_fliptorsion_2/')
#run_dir = Path('./testnewcode_2/')

# Ligand information
#ligand_reference = input_structures_dir / 'KEMP1_TSA_h.sdf'


Sculpt( 
        # folding
        folder=folder,
        fixers=fixers,

        # ga operators
        selector=selector,
        mutator=mutator,
        crossover=crossover,
        initial_population_mutator=initial_population_mutator,

        # scoring
        scoring_function=scoring_function,

        # params
        num_cycles=15,
        pop_size=30,
        
        # inputs
        structure_input_dir=input_structures_dir,
        run_dir=run_dir,
        sdf_files=[ligand_reference],
)

