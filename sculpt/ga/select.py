"""
Selection strategies for genetic algorithms using Design.score attribute.

All strategies accept a population of Design objects (each expected to have a
`.score` numeric attribute or None) and return a list of selected Designs of
length `k` (when possible).

This file defines a structural Protocol `SelectionStrategy` and several concrete
implementations: TopKSelection, TournamentSelection, ElitismSelection and
RouletteWheelSelection.
"""

from typing import Sequence, Protocol, List, TypeVar, Optional
import random
import numpy as np
from matplotlib import pyplot as plt

Design = TypeVar('Design')


class SelectionStrategy(Protocol[Design]):
    """Abstract selection strategy for genetic algorithms.

    Implementations should return a sequence of selected individuals from the
    provided population of length `k`.
    """

    def select(self, population: Sequence[Design], k: int) -> List[Design]:
        pass


def _score_value(d, highest: bool) -> float:
    """Helper to extract a finite numeric score from a design object.

    None scores are pushed to the worst end for deterministic behavior:
    - when highest=True, treat None as -inf so they lose to any numeric score
    - when highest=False, treat None as +inf so they lose when preferring low
    """
    s = getattr(d, 'score', None)
    if s is None:
        return float('-inf') if highest else float('inf')
    try:
        return float(s)
    except Exception:
        return float('-inf') if highest else float('inf')


def _prescale(fitness, fmultiple=2.0):
    """
    Compute linear scaling coefficients for genetic algorithm fitness scaling.

    Parameters
    ----------
    fitness : array-like
        Raw objective/fitness values.
    fmultiple : float
        Desired multiple of average fitness for the best individual.

    Returns
    -------
    a, b : float
        Coefficients for linear scaling: f_scaled = a * f + b
    """

    fitness = np.asarray(fitness)

    # remove nans
    fitness = fitness[~np.isnan(fitness)]

    if len(fitness) == 0:
        print('WARNING: All fitness values were NaN! This is a problem.')
        return 1.0, 0.0

    # if all fitness values are the same, return 1, 0
    if np.all(fitness == fitness[0]):
        print('WARNING: All fitness values are the same! This is a problem.')
        return 1.0, 0.0

    umax = fitness.max()
    umin = fitness.min()
    uavg = fitness.mean()

    print(umax, umin, uavg)

    # Non-negative test
    threshold = (fmultiple * uavg - umax) / (fmultiple - 1.0)

    if umin > threshold:
        # Normal scaling
        delta = umax - uavg
        a = (fmultiple - 1.0) * uavg / delta
        b = uavg * (umax - fmultiple * uavg) / delta
    else:
        # Scale as much as possible (ensure min fitness = 0)
        delta = uavg - umin
        a = uavg / delta
        b = -umin * uavg / delta

    return a, b


def _scale_fitness(fitness, a, b):
    """
    Apply the linear fitness scaling.

    Parameters
    ----------
    fitness : array-like
        Raw fitness values
    a, b : float
        Scaling coefficients

    Returns
    -------
    np.ndarray
        Scaled fitness values
    """
    fitness = np.asarray(fitness)
    scaled_fitness = a * fitness + b
    scaled_fitness[scaled_fitness < 0] = 0 # Because of floating point errors, we can sometimes dip below 0
    return scaled_fitness

# not finished

# class TopKSelection:
#     """Select the top-k individuals by score.

#     Args:
#         highest: If True (default) select highest scores; if False select lowest.
#     """

#     def __init__(self, highest: bool = True):
#         self.highest = highest

#     def select(self, population: Sequence[Design], k: int) -> List[Design]:
#         if k <= 0:
#             return []
#         reverse = self.highest
#         sorted_pop = sorted(
#             population,
#             key=lambda d: _score_value(d, self.highest),
#             reverse=reverse,
#         )

#         # Make copies:
#         returned_designs = [design.copy() for design in sorted_pop[:k]]
#         return returned_designs


# class TournamentSelection:
#     """Tournament selection using Design.score.

#     Args:
#         tournament_size: number of contenders per tournament (default 3).
#         with_replacement: whether to sample winners with replacement.
#         highest: True to prefer high scores, False to prefer low scores.
#     """

#     def __init__(self, tournament_size: int = 3, with_replacement: bool = True, highest: bool = True):
#         self.tournament_size = max(2, int(tournament_size))
#         self.with_replacement = with_replacement
#         self.highest = highest

#     def select(self, population: Sequence[Design], k: int) -> List[Design]:
#         if k <= 0:
#             return []
#         n = len(population)
#         if n == 0:
#             return []
#         selected: List[Design] = []
#         indices = list(range(n))
#         for _ in range(k):
#             contenders = random.sample(indices, min(self.tournament_size, len(indices)))
#             best_idx = max(contenders, key=lambda i: _score_value(population[i], self.highest))
#             selected.append(population[best_idx])
#             if not self.with_replacement:
#                 # remove the chosen index from the pool so it cannot be chosen again
#                 indices.remove(best_idx)
#                 if not indices:
#                     break

#         # Make copies:
#         returned_designs = [design.copy() for design in selected]
#         return returned_designs


# class ElitismSelection:
#     """Keep the top `elite_count` individuals (by score), then fill the rest
#     using a remainder strategy or random sampling.

#     Args:
#         elite_count: number of elites to keep.
#         remainder_strategy: optional SelectionStrategy used to fill the rest.
#         highest: True to prefer high scores when computing elites; False for low.
#     """

#     def __init__(self, elite_count: int, remainder_strategy: Optional[SelectionStrategy] = None, highest: bool = True):
#         self.elite_count = max(0, int(elite_count))
#         self.remainder_strategy = remainder_strategy
#         self.highest = highest

#     def select(self, population: Sequence[Design], k: int) -> List[Design]:
#         if k <= 0:
#             return []
#         n = len(population)
#         if n == 0:
#             return []
#         elite_k = min(self.elite_count, k)
#         paired = sorted(population, key=lambda d: _score_value(d, self.highest), reverse=self.highest)
#         elites = paired[:elite_k]
#         remaining = k - elite_k
#         if remaining <= 0:
#             # Make copies:
#             returned_designs = [design.copy() for design in elites]
#             return returned_designs
#         if self.remainder_strategy is None:
#             # default remainder: random sampling with replacement from non-elites
#             choices = paired if not elites else paired[elite_k:]
#             if not choices:
#                 # if there are no remaining choices (population smaller than elite_count)
#                 # Make copies:
#                 returned_designs = [design.copy() for design in elites]
#                 return returned_designs
            
#             # Make copies:
#             returned_designs = [design.copy() for design in elites] + [random.choice(choices) for _ in range(remaining)]
#             return returned_designs
#         else:
#             rest_pop = paired[elite_k:]
#             # If no rest_pop, allow remainder strategy to draw from full population
#             if not rest_pop:
#                 rest_pop = paired
#             # Make copies:
#             returned_designs = [design.copy() for design in elites] + self.remainder_strategy.select(rest_pop, remaining)
#             return returned_designs


class RouletteWheelSelection:
    """Fitness-proportionate selection using Design.score.

    Args:
        highest: True when larger scores are better (default). 
                If False, scores are treated as if they were 1/score. (This changes scaling!)
    """

    def __init__(self, highest: bool = True):
        self.highest = highest

    def select(self, population: Sequence[Design], k: int, plot_file: Optional[str] = None, return_reference: bool = False) -> List[Design]:
        """
        Selects k individuals from the population using roulette wheel selection.

        Args:   
            population: The population to select from.
            k: The number of individuals to select.
            plot_file: The file to save the pie chart to.
            return_reference: If True, returns the selected individuals by reference. If False, returns copies of the selected individuals.
        
        Returns:
            A list of k individuals selected with replacement from the population.
        """
        if k <= 0:
            return []
        n = len(population)
        if n == 0:
            return []
        # Build numeric scores; map None to 0 to avoid infinities for roulette.
        raw_scores = []
        for d in population:
            s = getattr(d, 'score', None)
            try:
                raw_scores.append(float(s) if s is not None else 0.0)
            except Exception:
                raw_scores.append(0.0)
        
        # If lower scores are better, invert the sign so selection prefers small values.
        # N.B. I expect this to act weird! Let's stick to positive fitness functions for now.
        if not self.highest:
            raw_scores = [1/s for s in raw_scores]

        scores = np.array(raw_scores)
        print('RAW SCORES:', scores)

        ### SCALING (Goldberg, pg. 79)
        a, b = _prescale(scores)
        print('COEFFICIENTS:', a, b)
        scaled_scores = _scale_fitness(scores, a, b)
        print('SCALED SCORES:', scaled_scores)
        total_score = np.sum(scaled_scores)
        
        # Are there nans?
        if np.isnan(scaled_scores).any():
            print("NaNs detected in scores")
            print('Normal scores:')
            print(scores)
            print('Scaled scores:')
            print(scaled_scores)
            # Filter out nans. For now, make them super low scoring:
            scaled_scores[np.isnan(scaled_scores)] = 0.0
            total_score = np.sum(scaled_scores)

            print('Scaled scores after filtering nans:')
            print(scaled_scores)

        # plot a pie chart of the scores:
        try:
            if plot_file is not None:
                plt.pie(scaled_scores, labels=[d.name for d in population])
                plt.savefig(plot_file)
                plt.close()
        except Exception as e:
            print(f"Error plotting pie chart: {e}")

        selected_indices = np.random.choice(len(population), size=k, p=scaled_scores/total_score)
        if return_reference:
            selected_designs = [population[i] for i in selected_indices]
        else:
            selected_designs = [population[i].copy() for i in selected_indices]

        return selected_designs
