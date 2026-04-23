from dataclasses import dataclass, field
from typing import List, Callable, Optional, Union
from pathlib import Path
import numpy as np
import dill

from sculpt.sculpt import Sculpt, _write_statistics
from sculpt.reporting.reporting import visualize_population_clusters

@dataclass
class SculptState:
    name: str 
    scoring_function: object 
    sdf_files: list = field(default_factory=list)
    weight: float = 1.0  
    fixers: list = field(default_factory=list)

def _write_multistate_statistics(generation_designs, parents_reference, run_dir: Path, cycle_dir: Path, cycle_idx: int, states: List[SculptState]):
    stats_csv = run_dir / 'multistate_run_statistics.csv'
    all_sequences_csv = run_dir / 'multistate_all_sequences.csv'
    
    write_header = not stats_csv.exists()
    
    with open(stats_csv, 'a') as f:
        if write_header:
            f.write('cycle,pop_size,avg_agg_fitness,')
            f.write(','.join([f"avg_{s.name},min_{s.name},max_{s.name}" for s in states]))
            f.write('\n')
            
        pop_size = len(generation_designs)
        agg_scores = [d.score for d in generation_designs if d.score is not None and not np.isnan(d.score)]
        avg_agg_fitness = np.mean(agg_scores) if agg_scores else float('nan')
        
        f.write(f'{cycle_idx},{pop_size},{avg_agg_fitness},')
        
        state_metrics = []
        for state in states:
            state_scores = [d.state_scores.get(state.name, float('nan')) for d in generation_designs if hasattr(d, 'state_scores')]
            clean_state_scores = [s for s in state_scores if s is not None and not np.isnan(s)]
            if clean_state_scores:
                state_metrics.append(f"{np.mean(clean_state_scores)},{np.min(clean_state_scores)},{np.max(clean_state_scores)}")
            else:
                state_metrics.append("nan,nan,nan")
                
        f.write(','.join(state_metrics) + '\n')

    write_header = not all_sequences_csv.exists()
    with open(all_sequences_csv, 'a') as f:
        if write_header:
            cols = 'cycle,design_name,parents,sequence,fitness,history,' + ','.join([f"score_{s.name}" for s in states])
            f.write(f'{cols}\n')
        for d in generation_designs:
            history_str = " / ".join([str(h) for h in d.history]).replace('\n', ' ')
            seq_val = "".join([line.strip() for line in d.sequence.splitlines() if not line.startswith(">")])
            parents_str = " / ".join([p.name for p in d.parents])
            
            st_scores = []
            if hasattr(d, 'state_scores'):
                st_scores = [str(d.state_scores.get(s.name, 'nan')) for s in states]
            else:
                st_scores = ['nan'] * len(states)
            
            f.write(f'{cycle_idx},{d.name},{parents_str},{seq_val},{d.score},{history_str},{",".join(st_scores)}\n')


class SculptMultistate(Sculpt):
    def __init__(
        self,
        folder = None,
        fixers: list = None,
        selector = None,
        mutator = None,
        crossover = None,
        initial_population_mutator = None,
        states: List[SculptState] = None,
        state_aggregator: Callable = None,
        num_cycles: int = None,
        pop_size: int = None,
        structure_input_dir = None,
        run_dir: str = None,
        resume_from_dir: str = None
    ):
        if resume_from_dir is not None:
            # The base handles loading and setting up everything via super() call below,
            # or we can directly invoke the logic. Wait, super() will handle resume_from_dir.
            pass

        if resume_from_dir is None:
            missing = []
            if folder is None: missing.append('folder')
            if selector is None: missing.append('selector')
            if mutator is None: missing.append('mutator')
            if crossover is None: missing.append('crossover')
            if not states: missing.append('states')
            if num_cycles is None: missing.append('num_cycles')
            if pop_size is None: missing.append('pop_size')
            if structure_input_dir is None: missing.append('structure_input_dir')
            if run_dir is None: missing.append('run_dir')
            
            if missing:
                raise ValueError(f"Missing required arguments for a new multistate run: {', '.join(missing)}")

            self.states = states
            # Default aggregator: Sum of weighted scores
            self.state_aggregator = state_aggregator or (lambda scores_dict: sum(st.weight * (scores_dict.get(st.name, np.nan) or np.nan) for st in self.states))
            
        # We pass dummy scoring_function to super to avoid its ValueError validation
        super().__init__(
            folder=folder,
            fixers=fixers,
            selector=selector,
            mutator=mutator,
            crossover=crossover,
            initial_population_mutator=initial_population_mutator,
            scoring_function="multistate_dummy",
            num_cycles=num_cycles,
            pop_size=pop_size,
            structure_input_dir=structure_input_dir,
            run_dir=run_dir,
            sdf_files=[],
            resume_from_dir=resume_from_dir
        )

    def run(self):
        """Main loop to run the Sculpt Multistate pipeline."""
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

            ### Step 0: Calculate Fitness
            for design in self.current_generation:
                design.state_scores = getattr(design, 'state_scores', {})
                design.state_structures = getattr(design, 'state_structures', {})

            for state in self.states:
                all_folds = self.folder.fold_batch(self.current_generation, ligand_sdf_files=state.sdf_files)
                
                for design, folds in zip(self.current_generation, all_folds):
                    applied_fixers = self.fixers + state.fixers
                    for fixer in applied_fixers:
                        folds = [fixer(fold) for fold in folds]
                    
                    if len(folds) == 0:
                        design.state_scores[state.name] = np.nan
                        continue
                    
                    for f in folds: 
                        f.score = state.scoring_function.score(f)

                    best_fold = max(folds, key=lambda f: f.score)
                    design.state_scores[state.name] = float(np.mean([f.score for f in folds]))
                    design.state_structures[state.name] = best_fold.structure
                    
                    design.structure = best_fold.structure
                    design.structure_file(cycle_dir / '0_Folded' / f'{design.name}_{state.name}.cif')

            # Aggregate scores over all states
            for design in self.current_generation:
                if any(np.isnan(val) for val in design.state_scores.values()):
                    design.score = np.nan
                else:
                    design.score = self.state_aggregator(design.state_scores)
                
                if self.states:
                    design.structure = design.state_structures.get(self.states[0].name)
                
                print(f"Aggregated score for {design.name}: {design.score:.3f} | Breakdown: {design.state_scores}")
                design.sequence_file(cycle_dir / '0_Folded' / f'{design.name}.fasta')

            ### Step 1: Selection / Reproduction
            parents_reference = self.selector.select(self.current_generation, len(self.current_generation), plot_file=cycle_dir / '1_Reproduction' / f'score_pie_chart_{i}.png')

            ### Step 2: Crossover
            next_generation = self.crossover.crossover(parents_reference)
            for design in next_generation:
                design.sequence_file(cycle_dir / '2_Crossover' / f'{design.name}.fasta')

            ### Step 3: Mutation 
            try:
                next_generation = self.mutator.mutate(next_generation, save_EM_structures_dir=cycle_dir / '3_Mutation', unique_naming=False)
            except TypeError:
                next_generation = self.mutator.mutate(next_generation, unique_naming=False)
            for design in next_generation:
                design.sequence_file(cycle_dir / '3_Mutation' / f'{design.name}.fasta')
            
            ### Step 4: Write Statistics & Plots
            try:
                _write_statistics(self.current_generation, parents_reference, self.run_dir, cycle_dir, i)
                _write_multistate_statistics(self.current_generation, parents_reference, self.run_dir, cycle_dir, i, self.states)

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
