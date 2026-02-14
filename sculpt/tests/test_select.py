from sculpt.select import SculptResidueSelector
from sculpt.geometry import Residue

class FakeStruct:
    def __init__(self, residues):
        self._res = residues
    def get_residues(self):
        return self._res


def test_basic_ops():
    r1 = Residue('A', 1)
    r2 = Residue('A', 2)
    r3 = Residue('B', 1)

    s1 = SculptResidueSelector([r1, r2])
    s2 = SculptResidueSelector([r2, r3])
    fs = FakeStruct([r1, r2, r3])

    assert s1(fs) == [r1, r2]
    assert s2(fs) == [r2, r3]
    assert (s1 | s2)(fs) == [r1, r2, r3]
    assert (s1 & s2)(fs) == [r2]
    assert (s1 - s2)(fs) == [r1]
    assert (s1 ^ s2)(fs) == [r1, r3]


def test_complement_and_all():
    r1 = Residue('A', 1)
    r2 = Residue('A', 2)
    r3 = Residue('B', 1)

    s1 = SculptResidueSelector([r1, r2])
    all_sel = SculptResidueSelector.from_all()
    fs = FakeStruct([r1, r2, r3])

    assert all_sel(fs) == [r1, r2, r3]
    assert (~s1)(fs) == [r3]
