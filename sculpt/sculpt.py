from sculpt.ga.select import TournamentSelection
from sculpt.tasks.optimize import SculptOptimizer 
from sculpt.tasks.resequence import SculptResequencer
from sculpt.tasks.fold import SculptFolder, SculptHydrogenAdder
from sculpt.tasks.score import SculptGeometricScoringFunction

from sculpt.ga.select import RouletteWheelSelection, TournamentSelection
from sculpt.ga.mutate import SimpleMutation, DirectedMutation
from sculpt.ga.crossover import SimpleCrossover, SpatialCrossover

from sculpt.reporting.reporting import visualize_population_clusters

from sculpt.id import generate_id
from sculpt.design import Design

from shim import StandardMolecule

from pathlib import Path
import random
import numpy as np
import dill

def _write_statistics(generation_designs, parents_reference, run_dir: Path, cycle_dir: Path, cycle_idx: int):
    stats_csv = run_dir / 'run_statistics.csv'
    all_sequences_csv = run_dir / 'all_sequences.csv'
    scores = [d.score for d in generation_designs if d.score is not None]
    pop_size = len(generation_designs)
    parent_scores = [d.score for d in parents_reference if d.score is not None]
    
    if scores:
        avg_fitness = np.nanmean(scores)
        min_fitness = np.nanmin(scores)
        max_fitness = np.nanmax(scores)
        std_fitness = np.nanstd(scores)
        sum_fitness = np.nansum(scores)

        # Calculate Selection Intensity:
        selection_intensity = (np.nanmean(parent_scores) - avg_fitness) / std_fitness
    else:
        avg_fitness = min_fitness = max_fitness = std_fitness = sum_fitness = float('nan')
        selection_intensity = float('nan')
        
    write_header = not stats_csv.exists()
    with open(stats_csv, 'a') as f:
        if write_header:
            f.write('cycle,pop_size,selection_intensity,avg_fitness,min_fitness,max_fitness,std_fitness,sum_fitness\n')
        f.write(f'{cycle_idx},{pop_size},{selection_intensity},{avg_fitness},{min_fitness},{max_fitness},{std_fitness},{sum_fitness}\n')

    # Write all sequences to a CSV file
    write_header = not all_sequences_csv.exists()
    with open(all_sequences_csv, 'a') as f:
        if write_header:
            f.write('cycle,design_name,parents,sequence,fitness,history\n')
        for d in generation_designs:
            # To help with CSV parsing, we're replacing commas with slashes in the history and parent names
            history_str = " / ".join([str(h) for h in d.history])
            history_str = history_str.replace('\n', ' ') # strip endlines:
            seq_val = "".join([line.strip() for line in d.sequence.splitlines() if not line.startswith(">")])
            parent_names = " / ".join([p.name for p in d.parents])
            f.write(f'{cycle_idx},{d.name},{parent_names},{seq_val},{d.score},{history_str}\n')
        
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

class Sculpt:
    def __init__(
        self,
        folder = None,
        fixers: list = None,
        selector = None,
        mutator = None,
        crossover = None,
        initial_population_mutator = None,
        scoring_function = None,
        num_cycles: int = None,
        pop_size: int = None,
        structure_input_dir = None,
        run_dir: str = None,
        sdf_files: list = None,
        resume_from_dir: str = None
    ):
        if resume_from_dir is not None:
            self._load_state(resume_from_dir)
            return

        missing = []
        if folder is None: missing.append('folder')
        if selector is None: missing.append('selector')
        if mutator is None: missing.append('mutator')
        if crossover is None: missing.append('crossover')
        if scoring_function is None: missing.append('scoring_function')
        if num_cycles is None: missing.append('num_cycles')
        if pop_size is None: missing.append('pop_size')
        if structure_input_dir is None: missing.append('structure_input_dir')
        if run_dir is None: missing.append('run_dir')
        
        if missing:
            raise ValueError(f"Missing required arguments for a new run: {', '.join(missing)}")

        self.folder = folder
        self.fixers = fixers if fixers is not None else []
        self.selector = selector
        self.mutator = mutator
        self.crossover = crossover
        self.initial_population_mutator = initial_population_mutator
        self.scoring_function = scoring_function
        self.num_cycles = num_cycles
        self.pop_size = pop_size
        self.structure_input_dir = Path(structure_input_dir) if structure_input_dir else None
        self.run_dir = Path(run_dir) if run_dir else None
        self.sdf_files = sdf_files if sdf_files is not None else []

        self.current_cycle = 0
        self.current_generation = []
        self.input_designs = []
        
        self.setup()

    def setup(self):
        """Dedicated function to set up all folders and initial population"""
        if self.structure_input_dir is not None:
            self.structure_input_dir = Path(self.structure_input_dir)
            
        input_designs = []
        if self.structure_input_dir and self.structure_input_dir.exists():
            for file in self.structure_input_dir.iterdir():
                if file.suffix in ['.pdb', '.cif', '.mmcif']:
                    input_design = Design(name=file.stem)
                    input_design.load_structure(file)
                    input_design.load_sequence_from_structure()
                    input_designs.append(input_design)
                if file.suffix in ['.fasta', '.fa']:
                    input_design = Design(name=file.stem)
                    input_design.load_sequence(file)
                    input_designs.append(input_design)

        copied_input_designs = []
        if self.pop_size and len(input_designs) > 0:
            for i in range(self.pop_size - len(input_designs)):
                new_input_design = random.choice(input_designs).copy(set_parent=True)
                copied_input_designs.append(new_input_design)
            
            copied_input_designs_with_mutations = []
            for design in copied_input_designs:
                if self.initial_population_mutator:
                    mutated_design = self.initial_population_mutator.mutate([design], unique_naming=True)[0]
                else: 
                    mutated_design = design
                copied_input_designs_with_mutations.append(mutated_design)

            input_designs.extend(copied_input_designs_with_mutations)
            
        print('Input structures:')
        for design in input_designs:
            print(' -', design.name)

        if self.run_dir is not None:
            self.run_dir.mkdir(parents=True, exist_ok=True)

        self.input_designs = input_designs
        self.current_generation = input_designs

        for design in self.input_designs:
            print(design.name)
            
        self._save_state()

    def run(self):
        """Main loop to run the Sculpt pipeline."""
        if self.run_dir is None:
            print("Cannot run without a designated run_dir")
            return
            
        for i in range(self.current_cycle, self.num_cycles):
            self.current_cycle = i
            cycle_dir = self.run_dir / f'cycle_{i}'
            cycle_dir.mkdir(parents=True, exist_ok=True)
            
            step_dirs = ['0_Folded', '1_Reproduction', '2_Crossover', '3_Mutation']
            for step in step_dirs:
                (cycle_dir / step).mkdir(parents=True, exist_ok=True)

            all_folds = self.folder.fold_batch(self.current_generation, ligand_sdf_files=self.sdf_files)
            
            for design, folds in zip(self.current_generation, all_folds):
                for fixer in self.fixers:
                    folds = [fixer(fold) for fold in folds]
                
                if len(folds) == 0:
                    design.score = np.nan
                    continue
                
                for f in folds: f.score = self.scoring_function.score(f)

                best_fold = max(folds, key=lambda f: f.score)
                design.score = np.mean([f.score for f in folds])
                design.structure = best_fold.structure
                print(f"Average score for {design.name}: {design.score}")

                design.structure_file(cycle_dir / '0_Folded' / f'{design.name}.cif')
                design.sequence_file(cycle_dir / '0_Folded' / f'{design.name}.fasta')

            parents_reference = self.selector.select(self.current_generation, len(self.current_generation), plot_file=cycle_dir / '1_Reproduction' / f'score_pie_chart_{i}.png')

            next_generation = self.crossover.crossover(parents_reference)
            for design in next_generation:
                design.sequence_file(cycle_dir / '2_Crossover' / f'{design.name}.fasta')

            try:
                next_generation = self.mutator.mutate(next_generation, save_EM_structures_dir=cycle_dir / '3_Mutation', unique_naming=False)
            except TypeError:
                next_generation = self.mutator.mutate(next_generation, unique_naming=False)
            for design in next_generation:
                design.sequence_file(cycle_dir / '3_Mutation' / f'{design.name}.fasta')
            
            try:
                _write_statistics(self.current_generation, parents_reference, self.run_dir, cycle_dir, i)

                fig, labels, metrics = visualize_population_clusters(
                    seqs=[d.sequence for d in self.current_generation],
                    fitness=[d.score for d in self.current_generation],
                    gen=i,
                    eps=0.06,
                    min_samples=3,
                    embed="pca",
                    outpath= cycle_dir / '1_Reproduction' / f'diversity_gen_{i}.png',
                )
            except Exception as e:
                print(f"Error writing statistics and plots: {e}")

            self.current_generation = next_generation
            self.current_cycle += 1
            self._save_state()

    def _save_state(self):
        if not self.run_dir:
            return
        state_file = self.run_dir / 'sculpt_state.pkl'
        with open(state_file, 'wb') as f:
            dill.dump(self.__dict__, f)

    def _load_state(self, run_dir):
        run_dir = Path(run_dir)
        state_file = run_dir / 'sculpt_state.pkl'
        if not state_file.exists():
            raise FileNotFoundError(f"State file not found at {state_file}. Cannot resume run.")
        with open(state_file, 'rb') as f:
            state = dill.load(f)
        self.__dict__.update(state)
        print(f"Resumed run from cycle {self.current_cycle}")
