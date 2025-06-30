### The following are the optimizer classes, used to optimize a structure file using EasyMD.
import tempfile
from pathlib import Path
import ribbon

class SculptOptimizer:
    """Optimizer class for structure optimization using custom restraints.
    
    This class provides functionality to optimize molecular structures using EasyMD
    with custom bonds, angles, and torsions as restraints.
    """
    
    def __init__(self, custom_bonds=[], custom_angles=[], custom_torsions=[]):
        """Initialize the SculptOptimizer with custom restraints.
        
        Args:
            custom_bonds: List of custom bond restraints to apply during optimization.
            custom_angles: List of custom angle restraints to apply during optimization.
            custom_torsions: List of custom torsion restraints to apply during optimization.
        """
        self.custom_bonds = custom_bonds
        self.custom_angles = custom_angles
        self.custom_torsions = custom_torsions

    def optimize(self, structure_file, output_file, sdf_files = []):
        """Use EasyMD to optimize the structure file using custom restraints.
        
        Args:
            structure_file: The input structure file to optimize.
            output_file: The output file to write the optimized structure to.
            sdf_files: A list of SDF files to use for ligand parameterization.
        """
        # Create a temp directory:
        temp_dir = tempfile.TemporaryDirectory()

        ribbon.EasyMD(
            input_file=file,
            output_prefix=str(Path(temp_dir.name) / "optimized"),
            ligand_files=sdf_files,
            duration=1,
            custom_bonds=self.custom_bonds,
            custom_torsions=self.custom_torsions,
            custom_angles=self.custom_angles,
            minimize_only=True,
        ).run()

        # Grab the output structure (prefix "_EM.pdb") and copy it to the output file:
        output_structure_file = str(Path(temp_dir.name) / "optimized_EM.pdb")
        with open(output_structure_file, 'r') as infile, open(output_file, 'w') as outfile:
            outfile.write(infile.read())
        # Clean up the temporary directory
        temp_dir.cleanup()
        
        return

    def __call__(self, structure_file, output_file, sdf_files = []):
        """Allow the optimizer to be called directly.
        
        Args:
            structure_file: The input structure file to optimize.
            output_file: The output file to write the optimized structure to.
            sdf_files: A list of SDF files to use for ligand parameterization.
            
        Returns:
            The result of the optimize method.
        """
        return self.optimize(structure_file, output_file, sdf_files) 



    
