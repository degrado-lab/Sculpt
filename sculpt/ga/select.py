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

Design = TypeVar('Design')


class SelectionStrategy(Protocol[Design]):
    """Abstract selection strategy for genetic algorithms.

    Implementations should return a sequence of selected individuals from the
    provided population of length `k`.
    """

    def select(self, population: Sequence[Design], k: int) -> List[Design]:
        ...


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


class TopKSelection:
    """Select the top-k individuals by score.

    Args:
        highest: If True (default) select highest scores; if False select lowest.
    """

    def __init__(self, highest: bool = True):
        self.highest = highest

    def select(self, population: Sequence[Design], k: int) -> List[Design]:
        if k <= 0:
            return []
        reverse = self.highest
        sorted_pop = sorted(
            population,
            key=lambda d: _score_value(d, self.highest),
            reverse=reverse,
        )
        return sorted_pop[:k]


class TournamentSelection:
    """Tournament selection using Design.score.

    Args:
        tournament_size: number of contenders per tournament (default 3).
        with_replacement: whether to sample winners with replacement.
        highest: True to prefer high scores, False to prefer low scores.
    """

    def __init__(self, tournament_size: int = 3, with_replacement: bool = True, highest: bool = True):
        self.tournament_size = max(2, int(tournament_size))
        self.with_replacement = with_replacement
        self.highest = highest

    def select(self, population: Sequence[Design], k: int) -> List[Design]:
        if k <= 0:
            return []
        n = len(population)
        if n == 0:
            return []
        selected: List[Design] = []
        indices = list(range(n))
        for _ in range(k):
            contenders = random.sample(indices, min(self.tournament_size, len(indices)))
            best_idx = max(contenders, key=lambda i: _score_value(population[i], self.highest))
            selected.append(population[best_idx])
            if not self.with_replacement:
                # remove the chosen index from the pool so it cannot be chosen again
                indices.remove(best_idx)
                if not indices:
                    break
        return selected


class ElitismSelection:
    """Keep the top `elite_count` individuals (by score), then fill the rest
    using a remainder strategy or random sampling.

    Args:
        elite_count: number of elites to keep.
        remainder_strategy: optional SelectionStrategy used to fill the rest.
        highest: True to prefer high scores when computing elites; False for low.
    """

    def __init__(self, elite_count: int, remainder_strategy: Optional[SelectionStrategy] = None, highest: bool = True):
        self.elite_count = max(0, int(elite_count))
        self.remainder_strategy = remainder_strategy
        self.highest = highest

    def select(self, population: Sequence[Design], k: int) -> List[Design]:
        if k <= 0:
            return []
        n = len(population)
        if n == 0:
            return []
        elite_k = min(self.elite_count, k)
        paired = sorted(population, key=lambda d: _score_value(d, self.highest), reverse=self.highest)
        elites = paired[:elite_k]
        remaining = k - elite_k
        if remaining <= 0:
            return elites
        if self.remainder_strategy is None:
            # default remainder: random sampling with replacement from non-elites
            choices = paired if not elites else paired[elite_k:]
            if not choices:
                # if there are no remaining choices (population smaller than elite_count)
                return elites
            return elites + [random.choice(choices) for _ in range(remaining)]
        else:
            rest_pop = paired[elite_k:]
            # If no rest_pop, allow remainder strategy to draw from full population
            if not rest_pop:
                rest_pop = paired
            return elites + self.remainder_strategy.select(rest_pop, remaining)


class RouletteWheelSelection:
    """Fitness-proportionate selection using Design.score.

    Args:
        highest: True when larger scores are better (default). If False, lower
                 scores are considered better and selection will prefer smaller
                 numeric values.
    """

    def __init__(self, highest: bool = True):
        self.highest = highest

    def select(self, population: Sequence[Design], k: int) -> List[Design]:
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
        if not self.highest:
            raw_scores = [-s for s in raw_scores]
        min_score = min(raw_scores)
        shifted = [s - min_score for s in raw_scores]
        total = sum(shifted)
        if total <= 0:
            # Fall back to uniform random sampling
            return [random.choice(list(population)) for _ in range(k)]
        cum = []
        c = 0.0
        for s in shifted:
            c += s
            cum.append(c)
        selected: List[Design] = []
        for _ in range(k):
            r = random.random() * total
            idx = next(i for i, v in enumerate(cum) if v > r)
            selected.append(population[idx])
        return selected
