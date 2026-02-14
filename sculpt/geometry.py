### The following are classes to define custom bonds, angles, and torsions
### With the str() function, they output a string representation which EasyMD accepts as input.
from dataclasses import dataclass

@dataclass
class Atom:
    chain: str
    residue: int
    name: str

    def __str__(self):
        return f"{self.chain}:{self.residue}:{self.name}"

@dataclass
class Residue:
    chain: str
    residue: int

    def __str__(self):
        return f"{self.chain}:{self.residue}"

@dataclass
class CustomBond:
    atom_1: Atom
    atom_2: Atom
    force_constant: float
    target_distance: float #Angstroms

    def __str__(self):
        return f"{self.atom_1},{self.atom_2},{self.force_constant},{self.target_distance}"

@dataclass
class CustomAngle:
    atom_1: Atom
    atom_2: Atom
    atom_3: Atom
    force_constant: float
    target_angle: float #Degrees

    def __str__(self):
        return f"{self.atom_1},{self.atom_2},{self.atom_3},{self.force_constant},{self.target_angle}"

@dataclass
class CustomTorsion:
    atom_1: Atom
    atom_2: Atom
    atom_3: Atom
    atom_4: Atom
    force_constant: float
    periodicity: int
    target_angle: float #Degrees

    def __str__(self):
        return (f"{self.atom_1},{self.atom_2},{self.atom_3},{self.atom_4},{self.force_constant},{self.periodicity},{self.target_angle}")
