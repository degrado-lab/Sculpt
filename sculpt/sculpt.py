from sculpt.tasks.optimize import SculptOptimizer 
from sculpt.tasks.resequence import SculptResequencer
from sculpt.tasks.fold import SculptFolder, SculptHydrogenAdder
from sculpt.tasks.score import SculptGeometricScoringFunction
from sculpt.id import generate_id
from sculpt.design import Design

from shim import StandardMolecule

from pathlib import Path
import random

def sculpt(optimizer: SculptOptimizer, 
            resequencer: SculptResequencer, 
            folder: SculptFolder, 
            scoring_function: SculptGeometricScoringFunction,
            structure_input_dir,
            num_cycles: int = 10,
            run_dir: str = None,
            sdf_files: list = [],
            top_designs_num: int = 5,
            cycle_retry: int = 5):
    """Main function to run the Sculpt pipeline for ligand design and optimization.
    This function orchestrates the complete Sculpt pipeline, including energy minimization,
    sequence design with LigandMPNN, structure prediction with CHAI-1, hydrogen addition,
    scoring, and iterative refinement over multiple cycles.

    Args:
        optimizer: Optimizer object for structure optimization.
        resequencer: Resequencer object for generating new sequences.
        folder: Folder object for structure prediction.
        scoring_function: Scoring function for evaluating designs.
    """

    # Grab the current input structures:
    structure_input_dir = Path(structure_input_dir)
    current_input_designs = []
    for file in structure_input_dir.iterdir():
        if file.suffix in ['.pdb', '.cif', '.mmcif']:
            input_design = Design(name=file.stem)
            input_design.load_structure(file)
            current_input_designs.append(input_design)
    
    # Create standardized input structures:
    standard_molecules = [StandardMolecule(structure_file=sdf_file) for sdf_file in sdf_files]
    hydrogens_on_input_files = False # Assume no hydrogens initially
    for design in current_input_designs:
        design.structure.standardize(standard_molecules=standard_molecules, use_hydrogens=hydrogens_on_input_files) # do we want True here?

    #Next, add hydrogens to the input structures:
    if not hydrogens_on_input_files:
        hydrogen_adder = SculptHydrogenAdder()
        for i, design in enumerate(current_input_designs):
            design_with_h = hydrogen_adder.add_hydrogens(design, sdf_files)
            current_input_designs[i] = design_with_h
    else:
        print('Input structures already have hydrogens. Skipping hydrogen addition.')
        
    # List input structures:
    print('Input structures:')
    for design in current_input_designs:
        print(' -', design.name)

    # Make the run directory if it doesn't exist:
    if run_dir is not None:
        run_dir.mkdir(parents=True, exist_ok=True)

    # Run data CSV:
    run_dir = Path(run_dir)
    run_data_csv = run_dir / 'run_data.csv'

    # Output the standardized structures to a new directory:
    standardized_dir = run_dir / 'standardized_inputs'
    standardized_dir.mkdir(parents=True, exist_ok=True)
    for design in current_input_designs:
        design.structure_file(standardized_dir / f'{design.name}_standardized.cif')

    # What's the score to beat?
    best_score = None

    for i in range(num_cycles):

        go_to_next_cycle = False

        # We'll keep on the same cycle, until we get a better score than the last best.
        # This makes sure we're spending energy on improving our best designs.

        print('-----------------------------------')
        print('---------- Cycle:', i, '-------------')
        print('-----------------------------------')

        cycle_data = {}

        # This will store all designs from each cycle in (Design, score) tuples:
        scored_designs = []

        # Make the directory for this cycle:
        cycle_dir = run_dir / f'cycle_{i}'
        cycle_dir.mkdir(parents=True, exist_ok=True)

        # Make the directories for each step: 
        step_dirs = ['1_Optimize', '2_Resequence', '3_Fold', '4_Score', '5_Select']
        for step in step_dirs:
            (cycle_dir / step).mkdir(parents=True, exist_ok=True)

        #1 Optimize Geometry
        optimized_designs = []
        for input_design in current_input_designs:
            # I don't know why we're getting errors here sometimes.
            try:
                optimized_design = optimizer.optimize(input_design, sdf_files) #? Creates a new design with name [input_design.name]_opt
                optimized_designs.append(optimized_design)
                #Output 
                optimized_design.structure_file(cycle_dir / '1_Optimize' / f'{optimized_design.name}.cif')
            except Exception as e:
                print(f'Error optimizing {input_design.name}: {e}')
                continue

        cycle_retry_count = 0
        selected_structures = []

        while not go_to_next_cycle:

            if cycle_retry_count >= cycle_retry:
                print('Reached maximum retries for this cycle. Moving to next cycle.')
                go_to_next_cycle = True
                if len(selected_structures) > 0: # IF we have to continue on, use our new (inferior) designs. This keeps us from getting in a loop.
                    current_input_designs = selected_structures
                break    

            #2 Resequence (We will rename the new designs)
            resequenced_designs = []
            for opt_design in optimized_designs:
                
                current_resequenced_designs = resequencer.resequence(opt_design, unique_naming=True) #These designs will have new names!
                resequenced_designs += current_resequenced_designs
                for reseq_design in current_resequenced_designs:
                    reseq_design.sequence_file(cycle_dir / '2_Resequence' / f'{reseq_design.name}.fasta')

            print(f'Generated {len(resequenced_designs)} resequenced designs.')
            
            #3 Fold (& add hydrogens)
            folded_design_groups = []
            for reseq_design in resequenced_designs:
                # Should this output the same object, or a copy?
                folded_design_group = folder.fold(reseq_design, ligand_sdf_files=sdf_files) #These designs will have new names!
                folded_design_groups.append(folded_design_group)
                for fold_design in folded_design_group:
                    fold_design.structure_file(cycle_dir / '3_Fold' / f'{fold_design.name}.cif')

            #4 Score
            scored_design_groups = [] # For the next cycle
            for folded_design_group in folded_design_groups:

                scored_design_group = []
                for folded_design in folded_design_group:
                    score = scoring_function.score(folded_design)
                    folded_design.score = score
                    scored_design_group.append(folded_design)

                scored_design_groups.append(scored_design_group)
            
            # Temporary. For each group of designs, take the average score and asssign it to all:
            subcycle_scored_designs = [] # New list for this subcycle
            for design_group in scored_design_groups:
                scores = []
                for design in design_group:
                    scores.append(design.score)
                avg_score = sum(scores) / len(scores)
                print(f'Average score for group {design_group[0].name} is {avg_score}')
                for design in design_group:
                    # Take the average or the individual score?

                    scored_designs.append( (design, avg_score) )
                    subcycle_scored_designs.append( (design, avg_score) )
                    #design_score_pairs.append( (design, design.score) )

            # 5 Select
            # For the next cycle, we will use the top X scoring structures from this cycle.
            scored_designs.sort(key=lambda x: x[1]) # Sort by score (lower is better)
            subcycle_scored_designs.sort(key=lambda x: x[1]) # Sort by score (lower is better)

            # Let's do an experiment where we select 3x the top designs, and then randomly pair down.
            # TODO: Do I want to keep this?
            selected_structures = [pair[0] for pair in scored_designs[:3*top_designs_num]] # Keep the top X structures for the next cycle
            random.shuffle(selected_structures)
            selected_structures = selected_structures[:top_designs_num] # Randomly select X

            print('Selected designs for next cycle:')
            for design in selected_structures:
                print(' -', design.name, 'Score:', design.score)
                design.structure_file(cycle_dir / '5_Select' / f'{design.name}.cif')

            # Write a CSV with info on each design:
            with open(run_data_csv, 'a') as f:
                if i == 0:
                    f.write('cycle,design_name,parent_design,score\n')
                for pair in subcycle_scored_designs:
                    design = pair[0]
                    score = pair[1]
                    parents_name = design.parents[0].parents[0].parents[0].name if design.parents else 'input'
                    f.write(f'{i},{design.name},{parents_name},{score}\n')

            # Are any of the designs better than the best we've seen?
            cycle_best_score = min([pair[1] for pair in subcycle_scored_designs])
            print('Best score this cycle:', cycle_best_score)
            if best_score is None or cycle_best_score < best_score:
                best_score = cycle_best_score
                print('New best score! Moving to next cycle.')
                current_input_designs = selected_structures
                go_to_next_cycle = True
            else:
                print('No improvement over best score of', best_score, '. Repeating cycle.')
                go_to_next_cycle = False

            cycle_retry_count += 1