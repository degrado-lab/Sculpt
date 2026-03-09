from sculpt.tasks.optimize import SculptOptimizer 
from sculpt.tasks.resequence import SculptResequencer
from sculpt.tasks.fold import SculptFolder, SculptHydrogenAdder
from sculpt.tasks.score import SculptGeometricScoringFunction

from sculpt.ga.select import RouletteWheelSelection
from sculpt.ga.mutate import SimpleMutation
from sculpt.ga.crossover import SimpleCrossover

from sculpt.reporting.reporting import visualize_population_clusters

from sculpt.id import generate_id
from sculpt.design import Design

from shim import StandardMolecule

from pathlib import Path
import random
import numpy as np

def _write_statistics(generation_designs, run_dir: Path, cycle_dir: Path, cycle_idx: int):
    stats_csv = run_dir / 'run_statistics.csv'
    all_sequences_csv = run_dir / 'all_sequences.csv'
    scores = [d.score for d in generation_designs if d.score is not None]
    pop_size = len(generation_designs)
    
    if scores:
        avg_fitness = np.mean(scores)
        min_fitness = np.min(scores)
        max_fitness = np.max(scores)
        std_fitness = np.std(scores)
        sum_fitness = np.sum(scores)
    else:
        avg_fitness = min_fitness = max_fitness = std_fitness = sum_fitness = float('nan')
        
    write_header = not stats_csv.exists()
    with open(stats_csv, 'a') as f:
        if write_header:
            f.write('cycle,pop_size,avg_fitness,min_fitness,max_fitness,std_fitness,sum_fitness\n')
        f.write(f'{cycle_idx},{pop_size},{avg_fitness},{min_fitness},{max_fitness},{std_fitness},{sum_fitness}\n')

    # Write all sequences to a CSV file
    write_header = not all_sequences_csv.exists()
    with open(all_sequences_csv, 'a') as f:
        if write_header:
            f.write('cycle,design_name,parents,sequence,fitness\n')
        for d in generation_designs:
            seq_val = "".join([line.strip() for line in d.sequence.splitlines() if not line.startswith(">")])
            f.write(f'{cycle_idx},{d.name},{[p.name for p in d.parents]},{seq_val},{d.score}\n')
        
    # Histogram of fitness values (raw data)
    with open(cycle_dir / 'fitness_values.csv', 'w') as f:
        f.write('design_name,parents,fitness\n')
        for d in generation_designs:
            f.write(f'{d.name},{[p.name for p in d.parents]},{d.score}\n')
            
    # Histogram of residue identity at each position
    clean_sequences = []
    for d in generation_designs:
        if d.sequence:
            seq_val = "".join([line.strip() for line in d.sequence.splitlines() if not line.startswith(">")])
            clean_sequences.append(seq_val)
            
    if clean_sequences:
        max_len = max(len(seq) for seq in clean_sequences)
        with open(cycle_dir / 'residue_frequencies.csv', 'w') as f:
            f.write('position,residue,count\n')
            for pos in range(max_len):
                res_at_pos = [seq[pos] for seq in clean_sequences if pos < len(seq)]
                counts = {}
                for res in res_at_pos:
                    counts[res] = counts.get(res, 0) + 1
                for res, count in counts.items():
                    f.write(f'{pos},{res},{count}\n')

def sculpt(optimizer: SculptOptimizer, 
            resequencer: SculptResequencer, 
            folder: SculptFolder, 
            scoring_function: SculptGeometricScoringFunction,
            structure_input_dir,
            num_cycles: int = 10,
            run_dir: str = None,
            sdf_files: list = [],
            top_designs_num: int = 5,
            pop_size: int = 20,
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
    input_designs = []
    for file in structure_input_dir.iterdir():
        if file.suffix in ['.pdb', '.cif', '.mmcif']:
            input_design = Design(name=file.stem)
            input_design.load_structure(file)
            input_design.load_sequence_from_structure()
            input_designs.append(input_design)
        if file.suffix in ['.fasta', '.fa']:
            input_design = Design(name=file.stem)
            input_design.load_sequence(file)
            input_designs.append(input_design)

    # We'll pad out the rest of the first generation with highly-mutated sequences, based on the input sequences:
    mutator = SimpleMutation(mutation_rate=1.0, mutations_per_sequence=30, temperature=2.0) # instead of 10, 1
    copied_input_designs = []
    for i in range(pop_size - len(input_designs)):
        # Need unique names
        new_input_design = random.choice(input_designs).copy()
        new_input_design.name = f'{new_input_design.name}_mut_{i}'
        copied_input_designs.append(new_input_design)
    
    copied_input_designs_with_mutations = []
    for design in copied_input_designs:
        mutated_design = mutator.mutate([design])[0]
        copied_input_designs_with_mutations.append(mutated_design)

    input_designs.extend(copied_input_designs_with_mutations)
    
    # # Create standardized input structures:
    # standard_molecules = [StandardMolecule(structure_file=sdf_file) for sdf_file in sdf_files]
    # hydrogens_on_input_files = False # Assume no hydrogens initially
    # for design in input_designs:
    #     design.structure.standardize(standard_molecules=standard_molecules, use_hydrogens=hydrogens_on_input_files) # do we want True here?

    # #Next, add hydrogens to the input structures:
    # if not hydrogens_on_input_files:
    #     hydrogen_adder = SculptHydrogenAdder()
    #     for i, design in enumerate(input_designs):
    #         design_with_h = hydrogen_adder.add_hydrogens(design, sdf_files)
    #         input_designs[i] = design_with_h
    # else:
    #     print('Input structures already have hydrogens. Skipping hydrogen addition.')


        
    # List input structures:
    print('Input structures:')
    for design in input_designs:
        print(' -', design.name)

    # Make the run directory if it doesn't exist:
    if run_dir is not None:
        run_dir.mkdir(parents=True, exist_ok=True)

    # Run data CSV:
    run_dir = Path(run_dir)
    run_data_csv = run_dir / 'run_data.csv'

    ### I think this is defunct now that I'm standardizing below
    # Output the standardized structures to a new directory:
    standardized_dir = run_dir / 'standardized_inputs'
    standardized_dir.mkdir(parents=True, exist_ok=True)
    for design in input_designs:
        design.sequence_file(standardized_dir / f'{design.name}.fasta')
        #design.structure_file(standardized_dir / f'{design.name}_standardized.cif')


    for design in input_designs:
        print(design.name)

    ### Loop is:
    # Start with population. 
    # 0. Calculate Fitness
    # 1. Selection / Reproduction
    # 2. Crossover
    # 3. Mutation (directed)
    # 4. Write Statistics
    
    current_generation = input_designs
    for i in range(num_cycles):
        # Make the directory for this cycle:
        cycle_dir = run_dir / f'cycle_{i}'
        cycle_dir.mkdir(parents=True, exist_ok=True)

        # Make the directories for each step: 
        step_dirs = ['0_Folded', '1_Reproduction', '2_Crossover', '3_Mutation']
        for step in step_dirs:
            (cycle_dir / step).mkdir(parents=True, exist_ok=True)

        ### Step 0: Calculate Fitness
        for design in current_generation:

            # Fold design, usually 5x structures per design
            folds = folder.fold(design, ligand_sdf_files=sdf_files)
            
            # Score each fold
            for f in folds: f.score = scoring_function.score(f)

            # Average the scores, but keep the best fold
            best_fold = max(folds, key=lambda f: f.score)
            design.score = np.mean([f.score for f in folds])
            design.structure = best_fold.structure
            print(f"Average score for {design.name}: {design.score}")

            # Write out folded structures:
            design.structure_file(cycle_dir / '0_Folded' / f'{design.name}.cif')
            design.sequence_file(cycle_dir / '0_Folded' / f'{design.name}.fasta')

        ### Step 1: Selection / Reproduction
        # TODO: Move the RouletteWheelSelection object to be passed in as an input
        selector = RouletteWheelSelection(highest=True)
        parents_reference = selector.select(current_generation, len(current_generation), plot_file=cycle_dir / '1_Reproduction' / f'score_pie_chart_{i}.png', return_reference=True) #returns a REFERENCE of the parents

        ### Step 2: Crossover
        # TODO: Move the SimpleCrossover object to be passed in as an input
        crossover = SimpleCrossover(crossover_rate=0.5)
        next_generation = crossover.crossover(parents_reference) # returns COPIES as children
        for design in next_generation:
            design.sequence_file(cycle_dir / '2_Crossover' / f'{design.name}.fasta')

        ### Step 3: Mutation 
        # TODO: Move the SimpleMutation object to be passed in as an input
        mutator = SimpleMutation(mutation_rate=0.5, mutations_per_sequence=1, temperature=1.0)
        next_generation = mutator.mutate(next_generation)        # modifies children IN PLACE
        for design in next_generation:
            design.sequence_file(cycle_dir / '3_Mutation' / f'{design.name}.fasta')
        
        ### Step 4: Write Statistics & Plots
        try:
            _write_statistics(current_generation, run_dir, cycle_dir, i)

            fig, labels, metrics = visualize_population_clusters(
                seqs=[d.sequence for d in current_generation],
                fitness=[d.score for d in current_generation],     # or None
                gen=i,
                eps=0.06,                       # tune this (see below)
                min_samples=3,
                embed="pca",                    # or "umap" if you install umap-learn
                outpath= cycle_dir / '1_Reproduction' / f'diversity_gen_{i}.png',
            )
        except Exception as e:
            print(f"Error writing statistics and plots: {e}")

        current_generation = next_generation