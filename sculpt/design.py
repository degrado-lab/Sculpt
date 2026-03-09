from pathlib import Path
from dataclasses import dataclass, field
from typing import Union, Optional
import tempfile
import io
import os
from contextlib import contextmanager

from shim import ShimStructure

"""
I think I can abstract this more. In my mind, each "design" is a sequence with potentially several folded structures. 
Or, should I have one sequence-structure pair? And group designs together somehow?

Regardless, I should ingest a sequence / structure from a file (like I do now), but immediately standardize it / check it. 
Then I can store it as a ShimStructure. I can also create a sparse ShimSequence to do similar things for FASTAs.
This would remove the headache I'm currently having with all the file handling etc. 
And I could store many annotated structures per Design if I desire.
"""
@dataclass
class Design:
    """Simplified Design object storing sequence and structure as text.
    We provide getters and setters to read/write using files or temp files.
    N.B. Currently, it doesn't check for PDB- or CIF-like structure formats."""

    name: str
    sequence: Union[str, None] = None      # Raw FASTA text
    #structure: Union[str, None] = None     # Raw PDB/CIF text
    structure: Union[ShimStructure, None] = None     # Structure Object
    score: Union[float, None] = None

    # Metadata:
    parents: list = field(default_factory=list)  # List of parent Design objects
    history: list = field(default_factory=list)  # List of steps taken to create this
    notes: Optional[str] = None  # Any additional notes

    # Additional helpers:
    #_structure_suffix: str = field(default='.cif', init=False)  # Default structure suffix

    def load_sequence(self, filepath: Union[str, Path]):
        """ Read sequence from a file and set the sequence attribute"""
        path = Path(filepath)
        if not path.exists():
            raise ValueError(f"File {filepath} does not exist.")
        self.sequence = path.read_text()

    def load_structure(self, filepath: Union[str, Path]):
        """ Read structure from a file and set the structure attribute"""
        self.structure = ShimStructure(structure_file=str(filepath), force=True)

    def load_sequence_from_structure(self):
        """ Populates sequence from a pre-loaded structure"""
        if self.structure is None:
            raise ValueError("No structure data available.")

        # TODO: Later, we should add a method to ShimStructure to get the sequence directly.
        # structure.to_fasta(outfile) will write a fasta file to outfile.
        # We can then read it back in:
        with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
            self.structure.to_fasta(tmp.name)
            self.sequence = Path(tmp.name).read_text()

    def sequence_file(self, filepath: Union[str, Path]):
        """ Write sequence to a file"""
        path = Path(filepath)
        if self.sequence is None:
            raise ValueError("No sequence data available.")
        path.write_text(self.sequence)

    def structure_file(self, filepath: Union[str, Path]):
        """ Write structure to a file"""
        path = Path(filepath)
        if self.structure is None:
            raise ValueError("No structure data available.")
        # Does the file have the right suffix? If not, warn and change it:
        if path.suffix.lower() == '.pdb':
            self.structure.to_pdb(str(path))
        elif path.suffix.lower() in ['.cif', '.mmcif']:
            self.structure.to_cif(str(path))
        else:
            raise ValueError("Structure file must have .pdb, .cif, or .mmcif suffix.")

    @contextmanager
    def temp_sequence_file(self, suffix=".fasta"):
        """Context manager: create a temp sequence file that is auto-deleted."""
        if self.sequence is None:
            raise ValueError("No sequence data available.")
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)  # we just want the name; we'll write separately
        path = Path(path)
        try:
            self.sequence_file(path)  # write out the sequence
            yield path
        finally:
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    @contextmanager
    def temp_structure_file(self, suffix=".cif"):
        """Context manager: create a temp structure file that is auto-deleted."""
        if self.structure is None:
            raise ValueError("No structure data available.")
        fd, path = tempfile.mkstemp(suffix=suffix)
        os.close(fd)  # we just want the name; we'll write separately
        path = Path(path)
        try:
            self.structure_file(path)  # write out the structure
            yield path
        finally:
            try:
                path.unlink()
            except FileNotFoundError:
                pass

    def copy(self, set_parent: bool = False):
        """ Create a copy of this Design object"""
        copy = Design(
            name=self.name,
            sequence=self.sequence,
            structure=self.structure,
            score=self.score,
            parents=[parent for parent in self.parents],
            history=self.history.copy(),
            notes=self.notes
        )
        if set_parent:
            copy.parents = [self]
        return copy

### Example usage:

# design = Design(name="example_design")

# design.load_sequence("example.fasta")
# design.load_structure("example.pdb")

# print(design.sequence)
# print(design.structure)

# design.sequence_file("output.fasta")
# design.structure_file("output.pdb")

# design.score = 9.5
# print(design.score)