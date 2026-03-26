from sculpt.design import Design
from sculpt.geometry import Atom
import shim.fix
from pathlib import Path

class SculptResidueFlipper:
    """A task to flip a symmetrical residue (like ASP) such that a specific atom (flip_atom) 
    is closest to a target atom (target_atom).
    
    This is useful for ensuring consistent orientation of symmetrical residues during optimization 
    or scoring.
    """

    def __init__(self, target_atom: Atom, flip_atom: Atom):
        """Initialize the SculptResidueFlipper.

        Args:
            target_atom: The Atom object whose proximity we're using as a target.
            flip_atom: The Atom object on the residue we want to ensure is closest to the target.
        """
        self.target_atom = target_atom
        self.flip_atom = flip_atom

    def flip(self, design: Design) -> Design:
        """Flip the residue in the given design structure.

        Args:
            design: The Design object to modify.

        Returns:
            A new Design object with the flipped residue.
        """
        # Create a new Design object (copy of the original)
        flipped_design = design.copy(set_parent=True)
        
        # Use a temporary structure file for the flip operation
        with design.temp_structure_file() as structure_file:
            # We'll use a temporary file for the output of the shim function
            # Since shim.fix.rename_atom_by_proximity takes file paths:
            import tempfile
            with tempfile.NamedTemporaryFile(suffix=Path(structure_file).suffix, delete=False) as outfile:
                outfile_path = outfile.name

            try:
                # Call the shim function
                # def rename_atom_by_proximity(infile: str, outfile: str, 
                #         target_chain: str, 
                #         target_resid: int,
                #         target_atom_name: str,
                #         renaming_chain: str,
                #         renaming_resid: int,
                #         closest_atom_name: str,
                #         )
                try:
                    shim.fix.rename_atom_by_proximity(
                        infile=str(structure_file),
                        outfile=str(outfile_path),
                        target_chain=self.target_atom.chain,
                        target_resid=self.target_atom.residue,
                        target_atom_name=self.target_atom.name,
                        renaming_chain=self.flip_atom.chain,
                        renaming_resid=self.flip_atom.residue,
                        closest_atom_name=self.flip_atom.name
                    )
                    
                    # Load the flipped structure back into the new design
                    flipped_design.load_structure(outfile_path)
                    
                    # Update history
                    flipped_design.history.append(f"Flipped residue {self.flip_atom.chain}:{self.flip_atom.residue} "
                                                f"to make {self.flip_atom.name} closer to "
                                                f"{self.target_atom.chain}:{self.target_atom.residue}:{self.target_atom.name}.")
                except Exception as e:
                    print(f"Error flipping residue {self.flip_atom.chain}:{self.flip_atom.residue}: {e}")
                    flipped_design = design
            finally:
                # Cleanup the temporary output file
                if Path(outfile_path).exists():
                    Path(outfile_path).unlink()

        return flipped_design

    def __call__(self, design: Design) -> Design:
        """Allow the flipper to be called directly.
        
        Args:
            design: The Design object to flip.
        """
        return self.flip(design)
