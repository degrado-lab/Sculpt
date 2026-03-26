"""
Mutation strategies for genetic algorithms using Design objects.

This file defines a structural Protocol `MutationStrategy` and a concrete
implementation: SimpleMutation, which uses BLOSUM62 matrices for substitution probabilities.
"""

import random
import numpy as np
from typing import Sequence, Protocol, List, TypeVar, Tuple

from sculpt.id import generate_id
from sculpt.tasks.fold import SculptFolder, SculptHydrogenAdder
from sculpt.tasks.optimize import SculptOptimizer
from sculpt.tasks.resequence import SculptResequencer

from shim import StandardMolecule
from pathlib import Path

Design = TypeVar('Design')


class MutationStrategy(Protocol[Design]):
    """Abstract mutation strategy for genetic algorithms.

    Implementations should return a sequence of mutated individuals from the
    provided population.
    """

    def mutate(self, population: Sequence[Design], **kwargs) -> List[Design]:
        pass


def _parse_fasta(fasta_text: str) -> Tuple[str, str]:
    """Helper to parse a simple FASTA string into header and sequence."""
    if not fasta_text:
        return "", ""
    lines = fasta_text.strip().split('\n')
    if not lines:
        return "", ""
    header = lines[0] if lines[0].startswith('>') else ">"
    seq = "".join([line.strip() for line in lines if not line.startswith('>')])
    if not seq and not lines[0].startswith('>'):
        seq = lines[0].strip()
    return header, seq


def _create_fasta(header: str, seq: str, line_length: int = 80) -> str:
    """Helper to create a FASTA string."""
    header = header if header.startswith('>') else f">{header}"
    seq_lines = [seq[i:i+line_length] for i in range(0, len(seq), line_length)]
    return f"{header}\n" + "\n".join(seq_lines)


AA_single_letter_codes = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
blosum62_array = np.array([[4.0, -1.0, -2.0, -2.0, 0.0, -1.0, -1.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0, 0.0, -3.0, -2.0, 0.0], 
                        [-1.0, 5.0, 0.0, -2.0, -3.0, 1.0, 0.0, -2.0, 0.0, -3.0, -2.0, 2.0, -1.0, -3.0, -2.0, -1.0, -1.0, -3.0, -2.0, -3.0], 
                        [-2.0, 0.0, 6.0, 1.0, -3.0, 0.0, 0.0, 0.0, 1.0, -3.0, -3.0, 0.0, -2.0, -3.0, -2.0, 1.0, 0.0, -4.0, -2.0, -3.0], 
                        [-2.0, -2.0, 1.0, 6.0, -3.0, 0.0, 2.0, -1.0, -1.0, -3.0, -4.0, -1.0, -3.0, -3.0, -1.0, 0.0, -1.0, -4.0, -3.0, -3.0], 
                        [0.0, -3.0, -3.0, -3.0, 9.0, -3.0, -4.0, -3.0, -3.0, -1.0, -1.0, -3.0, -1.0, -2.0, -3.0, -1.0, -1.0, -2.0, -2.0, -1.0], 
                        [-1.0, 1.0, 0.0, 0.0, -3.0, 5.0, 2.0, -2.0, 0.0, -3.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0, -1.0, -2.0, -1.0, -2.0], 
                        [-1.0, 0.0, 0.0, 2.0, -4.0, 2.0, 5.0, -2.0, 0.0, -3.0, -3.0, 1.0, -2.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0], 
                        [0.0, -2.0, 0.0, -1.0, -3.0, -2.0, -2.0, 6.0, -2.0, -4.0, -4.0, -2.0, -3.0, -3.0, -2.0, 0.0, -2.0, -2.0, -3.0, -3.0], 
                        [-2.0, 0.0, 1.0, -1.0, -3.0, 0.0, 0.0, -2.0, 8.0, -3.0, -3.0, -1.0, -2.0, -1.0, -2.0, -1.0, -2.0, -2.0, 2.0, -3.0], 
                        [-1.0, -3.0, -3.0, -3.0, -1.0, -3.0, -3.0, -4.0, -3.0, 4.0, 2.0, -3.0, 1.0, 0.0, -3.0, -2.0, -1.0, -3.0, -1.0, 3.0], 
                        [-1.0, -2.0, -3.0, -4.0, -1.0, -2.0, -3.0, -4.0, -3.0, 2.0, 4.0, -2.0, 2.0, 0.0, -3.0, -2.0, -1.0, -2.0, -1.0, 1.0], 
                        [-1.0, 2.0, 0.0, -1.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -2.0, 5.0, -1.0, -3.0, -1.0, 0.0, -1.0, -3.0, -2.0, -2.0], 
                        [-1.0, -1.0, -2.0, -3.0, -1.0, 0.0, -2.0, -3.0, -2.0, 1.0, 2.0, -1.0, 5.0, 0.0, -2.0, -1.0, -1.0, -1.0, -1.0, 1.0], 
                        [-2.0, -3.0, -3.0, -3.0, -2.0, -3.0, -3.0, -3.0, -1.0, 0.0, 0.0, -3.0, 0.0, 6.0, -4.0, -2.0, -2.0, 1.0, 3.0, -1.0], 
                        [-1.0, -2.0, -2.0, -1.0, -3.0, -1.0, -1.0, -2.0, -2.0, -3.0, -3.0, -1.0, -2.0, -4.0, 7.0, -1.0, -1.0, -4.0, -3.0, -2.0], 
                        [1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -2.0, 0.0, -1.0, -2.0, -1.0, 4.0, 1.0, -3.0, -2.0, -2.0], 
                        [0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0, -2.0, -1.0, 1.0, 5.0, -2.0, -2.0, 0.0], 
                        [-3.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -2.0, -2.0, -3.0, -2.0, -3.0, -1.0, 1.0, -4.0, -3.0, -2.0, 11.0, 2.0, -3.0], 
                        [-2.0, -2.0, -2.0, -3.0, -2.0, -1.0, -2.0, -3.0, 2.0, -1.0, -1.0, -2.0, -1.0, 3.0, -3.0, -2.0, -2.0, 2.0, 7.0, -1.0], 
                        [0.0, -3.0, -3.0, -3.0, -1.0, -2.0, -2.0, -3.0, -3.0, 3.0, 1.0, -2.0, 1.0, -1.0, -2.0, -2.0, 0.0, -3.0, -1.0, 4.0]])
uniform_array = np.ones_like(blosum62_array) / len(AA_single_letter_codes)


class SimpleMutation:
    """Sequence mutation using BLOSUM62 frequencies for a population of Design objects.

    Args:
        mutation_rate: Probability of a sequence undergoing mutation.
        mutations_per_sequence: Number of point mutations to introduce if selected.
        temperature: Controls how strictly to follow BLOSUM scores vs uniform random.
                     T=1 favors likely BLOSUM mutations; higher T becomes more uniform.
        fixed_residues: A string specifying fixed residues (default None). (e.g. "A1 A2 A3")
    """
    #TODO: All code here assumes 1 chain

    def __init__(self, mutation_rate: float = 0.5, mutations_per_sequence: int = 1, temperature: float = 1.0, fixed_residues: str = None):
        self.mutation_rate = max(0.0, min(1.0, float(mutation_rate)))
        self.mutations_per_sequence = max(1, int(mutations_per_sequence))
        self.temperature = max(0.01, float(temperature))
        self.fixed_residues = fixed_residues
        
        self.aa_list = AA_single_letter_codes
        self.aa_to_idx = {aa: i for i, aa in enumerate(self.aa_list)}
        
        # Precompute transition probabilities: P(new=j | old=i)
        # We use a softmax over the BLOSUM scores, carefully excluding self-mutation.
        self.transition_probs = np.zeros_like(blosum62_array, dtype=float)
        
        for i in range(len(self.aa_list)):
            # Scale scores by temperature
            scores = blosum62_array[i] / self.temperature
            
            # Mask out self-mutation
            scores_safe = np.copy(scores)
            scores_safe[i] = -np.inf # effectively zeroes its probability
            
            # Max score for numerical stability during exponentiation
            max_score = np.max(scores_safe)
            exp_scores = np.exp(scores_safe - max_score)
            
            # Ensure self-mutation remains strictly zero, then normalize
            exp_scores[i] = 0.0 
            self.transition_probs[i] = exp_scores / np.sum(exp_scores)

        ### DEBUG
        #print(self.transition_probs)

    def mutate(self, population: Sequence[Design], unique_naming: bool = False) -> List[Design]:
        """Perform mutation on a population.
        
        Args:
            population: Sequence of Design objects to mutate.
            unique_naming: If True, generate unique names for the offspring.
                           Otherwise, derive names from parents.
        
        Returns:
            A new list of mutated Design objects.
        """
        mutated_pop = []
        for parent in population:
            # We copy each parent directly, assuming they are either from 
            # parents_reference or the offspring copies generated by crossover.
            child = parent.copy() # Making exact copies, keeping original parents
            
            # Decide if this child mutates
            mutated = False
            mutations_made = []
            
            if random.random() < self.mutation_rate and child.sequence:
                h, seq = _parse_fasta(child.sequence)
                
                if len(seq) > 0:
                    seq_chars = list(seq)
                    
                    # Choose which positions to mutate
                    num_muts = min(self.mutations_per_sequence, len(seq_chars))
                    fixed_indices = {int("".join(filter(str.isdigit, r))) - 1 for r in self.fixed_residues.split()} if self.fixed_residues else set()
                    eligible_indices = [i for i in range(len(seq_chars)) if i not in fixed_indices]
                    num_muts = min(self.mutations_per_sequence, len(eligible_indices))
                    pos_to_mutate = random.sample(eligible_indices, num_muts)
                    
                    for pos in pos_to_mutate:
                        orig_aa = seq_chars[pos]
                        
                        # If it's a standard amino acid, use our parameterized BLOSUM probabilities
                        if orig_aa in self.aa_to_idx:
                            orig_idx = self.aa_to_idx[orig_aa]
                            probs = self.transition_probs[orig_idx]
                            new_aa = np.random.choice(self.aa_list, p=probs)
                        else:
                            # Fallback if non-standard: uniform over all 20 std amino acids
                            new_aa = np.random.choice(self.aa_list)

                        ### DEBUG
                        #print('MUTATING')
                        #print(orig_aa, new_aa)  
                        
                        # Only track if an actual change was made (though probabilities exclude self-mutation anyway)
                        if orig_aa != new_aa:
                            seq_chars[pos] = new_aa
                            mutations_made.append(f"{orig_aa}{pos+1}{new_aa}")
                    
                    if len(mutations_made) > 0:
                        mutated = True
                        new_seq = "".join(seq_chars)
                        
                        # If we don't use unique_naming, we just keep the same name as the parent.
                        # This is because we rename in the "crossover" step. Mutating is just a 
                        # modification of the parent sequence, not a new design.
                        if unique_naming:
                            c_name = generate_id()
                        else:
                            c_name = f"{parent.name}"
                            
                        child.name = c_name
                        child.sequence = _create_fasta(c_name, new_seq)
                        child.structure = None
                        
                        muts_str = ", ".join(mutations_made)
                        child.history.append(f"Mutation(s): {muts_str}")
                    
            mutated_pop.append(child)
            
        # Sanity check
        for child in mutated_pop:
            print(child.sequence)
        return mutated_pop

    def __call__(self, population: Sequence[Design], unique_naming: bool = False) -> List[Design]:
        """Allow the mutation strategy to be called directly.
        
        Args:
            population: Sequence of Design objects to mutate.
            unique_naming: If True, generate unique names for the offspring.
            
        Returns:
            The result of the mutate method.
        """
        return self.mutate(population, unique_naming=unique_naming)

class DirectedMutation:
    """Sequence mutation using user-provided force constraints in energy minimization, followed by MPNN.
    Structure is folded first, if necessary.

    Args:
        mutation_rate: Probability of a sequence undergoing mutation.
        resequencer: Resequencer object for generating new sequences.
        optimizer: Optimizer object for structure optimization.
        folder: Folder object for structure folding.
        sdf_files: List of SDF files for ligand parameterization.
        refold: If True, refold the structure even if it exists.
    """

    def __init__(
        self,
        mutation_rate: float = 0.5,
        resequencer: SculptResequencer = None,
        optimizer: SculptOptimizer = None,
        folder: SculptFolder = None,
        hydrogen_adder: SculptHydrogenAdder = None,
        refold: bool = False,
        sdf_files: List[str] = None,
    ):
        self.mutation_rate = mutation_rate
        self.resequencer = resequencer
        self.optimizer = optimizer
        self.folder = folder
        self.hydrogen_adder = hydrogen_adder
        self.refold = refold
        self.sdf_files = sdf_files

    def mutate(self, population: Sequence[Design], unique_naming: bool = False, save_EM_structures_dir: str = None) -> List[Design]:
        """Perform directed mutation on a population.
        
        Args:
            population: Sequence of Design objects to mutate.
            unique_naming: If True, generate unique names for the offspring.
            save_EM_structures_dir: Directory to save the EM structures of the mutated designs.
        Returns:
            A new list of mutated (or original) Design objects.
        """
        

        mutated_pop = []
        for parent in population:
            # Decide if this child mutates
            if random.random() < self.mutation_rate:
                # Start with a copy
                child = parent.copy()
                
                # 1) Fold/Hydrogens
                if self.refold or child.structure is None:
                    folder = self.folder
                    folded_results = folder.fold(child, ligand_sdf_files=self.sdf_files)
                    if folded_results:
                        child = folded_results[0]
                elif self.hydrogen_adder:
                    h_adder = self.hydrogen_adder
                    child = h_adder.add_hydrogens(child, ligand_sdf_files=self.sdf_files)

                # DEBUG
                # # write temporary structure file:
                # child.structure_file("temp_structure.pdb")
                # standard_molecules = [StandardMolecule(structure_file=sdf_file) for sdf_file in self.sdf_files]
                # hydrogens_on_input_files = False # Assume no hydrogens initially
                # child.structure.standardize(standard_molecules=standard_molecules, use_hydrogens=hydrogens_on_input_files)

                # 2) Energy-minimize with user-provided constraints
                optimizer = self.optimizer
                child = optimizer.optimize(child, sdf_files=self.sdf_files, unique_naming=unique_naming)
                
                if save_EM_structures_dir:
                    save_path = Path(save_EM_structures_dir) / f"{child.name}_EM.pdb"
                    child.structure_file(save_path)

                # 3) Run LigandMPNN/LaserMPNN to generate the new sequence
                resequencer = self.resequencer
                # Resequencer returns a list of designs
                reseq_results = resequencer.resequence(child, unique_naming=unique_naming)
                if reseq_results:
                    child = reseq_results[0]
                    child.history.append(f"DirectedMutation applied (model={self.folder.model}, rate={self.mutation_rate})")
                
                mutated_pop.append(child)
            else:
                # No mutation, add a copy of the parent to the new population
                mutated_pop.append(parent.copy())
        
        return mutated_pop

    def __call__(self, population: Sequence[Design], unique_naming: bool = False,  save_EM_structures_dir: str = None) -> List[Design]:
        """Allow the mutation strategy to be called directly.
        
        Args:
            population: Sequence of Design objects to mutate.
            unique_naming: If True, generate unique names for the offspring.
            
        Returns:
            The result of the mutate method.
        """
        return self.mutate(population, unique_naming=unique_naming)


class SimpleAndDirectedMutation:
    """Sequence mutation using both Simple and Directed mutation strategies.
    
    Applies simple mutation first, followed by directed mutation, using their
    respective mutation rates.

    Args:
        simple_mutation_rate: Probability of a sequence undergoing simple mutation.
        directed_mutation_rate: Probability of a sequence undergoing directed mutation.
        mutations_per_sequence: Number of point mutations to introduce if selected for simple mutation.
        temperature: Controls how strictly to follow BLOSUM scores vs uniform random for simple mutation.
        fixed_residues: A string specifying fixed residues for simple mutation (default None).
        resequencer: Resequencer object for generating new sequences in directed mutation.
        optimizer: Optimizer object for structure optimization in directed mutation.
        folder: Folder object for structure folding in directed mutation.
        hydrogen_adder: HydrogenAdder object for directed mutation.
        refold: If True, refold the structure even if it exists in directed mutation.
        sdf_files: List of SDF files for ligand parameterization in directed mutation.
    """

    def __init__(
        self,
        simple_mutation_rate: float = 0.5,
        directed_mutation_rate: float = 0.5,
        mutations_per_sequence: int = 1,
        temperature: float = 1.0,
        fixed_residues: str = None,
        resequencer: SculptResequencer = None,
        optimizer: SculptOptimizer = None,
        folder: SculptFolder = None,
        hydrogen_adder: SculptHydrogenAdder = None,
        refold: bool = False,
        sdf_files: List[str] = None,
    ):
        self.simple_mutation_rate = simple_mutation_rate
        self.directed_mutation_rate = directed_mutation_rate
        self.simple_mutation = SimpleMutation(
            mutation_rate=simple_mutation_rate,
            mutations_per_sequence=mutations_per_sequence,
            temperature=temperature,
            fixed_residues=fixed_residues
        )
        self.directed_mutation = DirectedMutation(
            mutation_rate=directed_mutation_rate,
            resequencer=resequencer,
            optimizer=optimizer,
            folder=folder,
            hydrogen_adder=hydrogen_adder,
            refold=refold,
            sdf_files=sdf_files
        )

    def mutate(self, population: Sequence[Design], unique_naming: bool = False, save_EM_structures_dir: str = None) -> List[Design]:
        """Perform both simple and directed mutation sequentially on a population.
        
        Args:
            population: Sequence of Design objects to mutate.
            unique_naming: If True, generate unique names for the offspring.
            save_EM_structures_dir: Directory to save the EM structures of the sequentially mutated designs.
            
        Returns:
            A new list of mutated (or original) Design objects.
        """
        # First apply simple mutation
        pop_after_simple = self.simple_mutation.mutate(population, unique_naming=unique_naming)
        
        # Then apply directed mutation
        final_pop = self.directed_mutation.mutate(pop_after_simple, unique_naming=unique_naming, save_EM_structures_dir=save_EM_structures_dir)
        
        return final_pop

    def __call__(self, population: Sequence[Design], unique_naming: bool = False, save_EM_structures_dir: str = None) -> List[Design]:
        """Allow the mutation strategy to be called directly.
        
        Args:
            population: Sequence of Design objects to mutate.
            unique_naming: If True, generate unique names for the offspring.
            save_EM_structures_dir: Directory to save the EM structures of the sequentially mutated designs.
            
        Returns:
            The result of the mutate method.
        """
        return self.mutate(population, unique_naming=unique_naming, save_EM_structures_dir=save_EM_structures_dir)
