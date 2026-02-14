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

from geometry import Residue

class SculptResidueSelector:
    """This class allows for creating on-the-fly residue selections based on a user-defined criteria.
    """
    
    def __init__(self, residues_list=[]):
        """Initialize the residue selector.
        
        Args:
            residues_list: A list of Residue objects to select.
        """
        self.residues_list = residues_list
        self.inverted = False               # Track if the selection is inverted. We need to hold on to this, until we're called with a structure file.

    def select(self, structure_file):
        """Select residues from a structure file.
        
        Args:
            structure_file: Path to the structure file to select residues from.
            
        Returns:
            A list of selected residues.
        """
        # If residues_list contains "all", return all residues.
        if 'all' in self.residues_list:
            # Use ribbon to get all residues from the structure file.
            import ribbon
            structure = ribbon.load_structure(structure_file)
            self.residues_list = [Residue(chain=res.chain, residue=res.residue) for res in structure.get_residues()]

        if self.inverted:
            all_residues_selector = self._get_all_residues(structure_file)
            self.residues_list = (all_residues_selector - self).residues_list

        # Call the selection function ("function") with the structure file
        return self.residues_list

    def __call__(self, structure_file):
        """Allow the selection function to be called directly.
        
        Args:
            structure_file: Path to the structure file to select residues from.
            
        Returns:
            A list of selected residues.
        """
        return self.select(structure_file)
      
    def _get_all_residues(self, structure_file):
        """Helper function to get all residues from a structure file.
        
        Args:
            structure_file: Path to the structure file. 

        Returns:
            A SculptResidueSelector object with a list of all Residues in the structure file.
        """
        # Need some method of opening the structure file and getting all residues.
        # Probably a job for Shim.
        # For now, just return a dummy list.
        return SculptResidueSelector(residues_list=[Residue(chain='A', residue=1), Residue(chain='A', residue=2), Residue(chain='B', residue=1), Residue(chain='B', residue=2)])

    ### ------------ LOGIC ------------ ###
    def __or__(self, other):
        """Combine two selectors with a logical OR.
        
        Args:
            other: Another SculptResidueSelector.
            
        Returns:
            A new SculptResidueSelector that selects residues from either selector.
        """
        if not isinstance(other, SculptResidueSelector):
            raise ValueError("Can only combine with another SculptResidueSelector.")
        
        combined_residues = list(set(self.residues_list + other.residues_list))
        return SculptResidueSelector(residues_list=combined_residues)

    def __and__(self, other):
        """Combine two selectors with a logical AND.
        
        Args:
            other: Another SculptResidueSelector.
            
        Returns:
            A new SculptResidueSelector that selects residues from both selectors.
        """
        if not isinstance(other, SculptResidueSelector):
            raise ValueError("Can only combine with another SculptResidueSelector.")
        
        common_residues = [res for res in self.residues_list if res in other.residues_list]
        return SculptResidueSelector(residues_list=common_residues)

    def __sub__(self, other):
        """Subtract residues from another selector.
        
        Args:
            other: Another SculptResidueSelector.
        """
        if not isinstance(other, SculptResidueSelector):
            raise ValueError("Can only combine with another SculptResidueSelector.")

        unique_residues = [res for res in self.residues_list if res not in other.residues_list]
        return SculptResidueSelector(residues_list=unique_residues)

    def __add__(self, other):
        """Add residues from another selector.
        
        Args:
            other: Another SculptResidueSelector.
        """
        return self.__or__(other)

    def __invert__(self):
        """Invert the selection (select all residues not in the current selection).
        
        Returns:
            A new SculptResidueSelector that selects all residues not in the current selection.
        """
        # This requires knowledge of all residues in the structure file.
        # For now, just raise an error.
        all_residues_selector = self._get_all_residues()
        return all_residues_selector - self


# Example usage:
if __name__ == "__main__":
    
    # Needs to work in like this:
    selector1 = SculptResidueSelector(residues_list=[Residue(chain='A', residue=1), Residue(chain='A', residue=2)])
    selector2 = SculptResidueSelector(residues_list=[Residue(chain='A', residue=2), Residue(chain='B', residue=1)])

    inverted_selector = ~selector1

    combined_selector = selector1 | inverted_selector

    combined_selector('test_file.pdb')
    