# WIP

"""
The idea here is to have classes that can on-the-fly select residues, atoms, etc. from a structure file based on a user-defined function.
So rather than manually defining residues to redesign with LigandMPNN, for example, we could define a function 
that selects residues based on proximity to the ligand, or based on some other criteria.

This would allow us to be flexible and do things like on-the-fly loop hallucinations using RFDiffusion etc. Because the overall # of residues wouldn't matter.
"""

"""
Somethings not quite right here. I have to think about the delayed logic more.
"""

from dataclasses import dataclass
from typing import Iterable, Callable, Optional, Set, Union
from geometry import Residue


@dataclass(frozen=True)
class _ResidueKey:
    chain: str
    residue: int

    @classmethod
    def from_residue(cls, r: Residue):
        return cls(chain=r.chain, residue=r.residue)


class SculptResidueSelector:
    """Lazy selector for residues supporting set-like logical operations.

    Instances are immutable: logical operators build expression trees that are
    evaluated only when the selector is called with a structure (path or object).

    Supported operators: | (union), & (intersection), - (difference), ^ (xor),
    ~ (complement). The complement operator produces an expression node that
    needs the full set of residues from the target structure to evaluate.

    A selector can be constructed from:
    - an iterable of `Residue` objects (concrete selection)
    - the special string 'all' to represent all residues (lazy until eval)
    - a callable that takes a structure and returns an iterable of `Residue`

    Example:
        s1 = SculptResidueSelector([Residue('A', 1), Residue('A', 2)])
        s2 = SculptResidueSelector([Residue('A', 2), Residue('B', 1)])
        combined = (s1 | s2) & ~SculptResidueSelector('all')
        selected = combined(structure)  # evaluate against structure
    """

    # Expression node types
    _OP_UNION = 'union'
    _OP_INTERSECT = 'intersect'
    _OP_DIFF = 'diff'
    _OP_XOR = 'xor'
    _OP_COMPLEMENT = 'complement'
    _OP_LEAF = 'leaf'

    def __init__(self,
                 residues: Optional[Union[Iterable[Residue], str, Callable]] = None,
                 _op: Optional[str] = None,
                 _children: Optional[Iterable['SculptResidueSelector']] = None):
        # Internal representation: an expression tree node
        if _op is None:
            # Leaf
            if residues is None:
                residues = []
            self._op = self._OP_LEAF
            if isinstance(residues, str) and residues == 'all':
                self._leaf_all = True
                self._leaf_set: Set[_ResidueKey] = set()
                self._leaf_func = None
            elif callable(residues):
                # user-provided function to compute residues when evaluating
                self._leaf_all = False
                self._leaf_set = set()
                self._leaf_func = residues
            else:
                self._leaf_all = False
                self._leaf_func = None
                self._leaf_set = set(_ResidueKey.from_residue(r) for r in (residues or []))
            self._children = None
        else:
            # Internal op node
            self._op = _op
            self._children = tuple(_children or ())
            self._leaf_all = False
            self._leaf_set = set()
            self._leaf_func = None

    # Construction helpers
    @classmethod
    def from_all(cls):
        return cls('all')

    @classmethod
    def from_callable(cls, func: Callable):
        return cls(func)

    # Evaluation
    def __call__(self, structure) -> list:
        """Evaluate the selector against a structure.

        `structure` can be a path-like object (string) or an already-loaded
        structure object that provides `get_residues()` returning Residue-like
        objects with `.chain` and `.residue` attributes.
        """
        keys = self._eval(structure)
        # Convert back to Residue dataclass instances. Return a list (Residue is
        # not hashable in this codebase) and sort for deterministic output.
        residues = [Residue(chain=k.chain, residue=k.residue) for k in keys]
        residues.sort(key=lambda r: (r.chain, r.residue))
        return residues

    def _eval(self, structure) -> Set[_ResidueKey]:
        # Leaf
        if self._op == self._OP_LEAF:
            if self._leaf_all:
                return self._get_all_keys(structure)
            if self._leaf_func is not None:
                residues = self._leaf_func(structure)
                return set(_ResidueKey.from_residue(r) for r in residues)
            return set(self._leaf_set)

        # Internal ops: evaluate children first
        if self._op == self._OP_UNION:
            out = set()
            for c in self._children:
                out |= c._eval(structure)
            return out
        if self._op == self._OP_INTERSECT:
            children = [c._eval(structure) for c in self._children]
            if not children:
                return set()
            out = children[0]
            for s in children[1:]:
                out &= s
            return out
        if self._op == self._OP_DIFF:
            a, b = self._children
            return a._eval(structure) - b._eval(structure)
        if self._op == self._OP_XOR:
            a, b = self._children
            return a._eval(structure) ^ b._eval(structure)
        if self._op == self._OP_COMPLEMENT:
            child = self._children[0]
            all_keys = self._get_all_keys(structure)
            return all_keys - child._eval(structure)
        raise RuntimeError(f"Unknown op {self._op}")

    def _get_all_keys(self, structure) -> Set[_ResidueKey]:
        # Accept either a path string or an already loaded object with get_residues()
        # If given a string, try to import ribbon (best-effort) and load it.
        residues = []
        if isinstance(structure, str):
            try:
                import ribbon
                st = ribbon.load_structure(structure)
                residues = list(st.get_residues())
            except Exception:
                # Fallback: empty
                residues = []
        else:
            # assume structure-like
            try:
                residues = list(structure.get_residues())
            except Exception:
                # Can't iterate; assume empty
                residues = []
        return set(_ResidueKey.from_residue(r) for r in residues)

    # Logical operators build new nodes
    def __or__(self, other: 'SculptResidueSelector') -> 'SculptResidueSelector':
        if not isinstance(other, SculptResidueSelector):
            return NotImplemented
        return SculptResidueSelector(_op=self._OP_UNION, _children=(self, other))

    def __and__(self, other: 'SculptResidueSelector') -> 'SculptResidueSelector':
        if not isinstance(other, SculptResidueSelector):
            return NotImplemented
        return SculptResidueSelector(_op=self._OP_INTERSECT, _children=(self, other))

    def __sub__(self, other: 'SculptResidueSelector') -> 'SculptResidueSelector':
        if not isinstance(other, SculptResidueSelector):
            return NotImplemented
        return SculptResidueSelector(_op=self._OP_DIFF, _children=(self, other))

    def __xor__(self, other: 'SculptResidueSelector') -> 'SculptResidueSelector':
        if not isinstance(other, SculptResidueSelector):
            return NotImplemented
        return SculptResidueSelector(_op=self._OP_XOR, _children=(self, other))

    def __invert__(self) -> 'SculptResidueSelector':
        return SculptResidueSelector(_op=self._OP_COMPLEMENT, _children=(self,))

    # Convenience for viewing/debugging
    def __repr__(self):
        if self._op == self._OP_LEAF:
            if self._leaf_all:
                return "SculptResidueSelector(all)"
            if self._leaf_func:
                return f"SculptResidueSelector(func={self._leaf_func})"
            return f"SculptResidueSelector({sorted(self._leaf_set, key=lambda k: (k.chain,k.residue))})"
        return f"SculptResidueSelector(op={self._op}, children={self._children})"


# Example usage (kept for developer reference)
if __name__ == "__main__":
    s1 = SculptResidueSelector([Residue(chain='A', residue=1), Residue(chain='A', residue=2)])
    s2 = SculptResidueSelector([Residue(chain='A', residue=2), Residue(chain='B', residue=1)])

    inverted = ~s1
    combined = s1 | inverted

    # The following will attempt to evaluate; if 'test.pdb' isn't loadable
    # by `ribbon`, the result will be empty.
    print(combined('test.pdb'))

    