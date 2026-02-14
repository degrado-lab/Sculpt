### The following are the optimizer classes, used to optimize a structure file using EasyMD.
from typing import List, Union

import tempfile
from pathlib import Path
import ribbon

from sculpt.design import Design
from shim import StandardMolecule

class SculptOptimizer:
    """Optimizer class for structure optimization using custom restraints.
    
    This class provides functionality to optimize molecular structures using EasyMD
    with custom bonds, angles, and torsions as restraints.
    """
    
    def __init__(self, custom_bonds=[], custom_angles=[], custom_torsions=[], full_sim: bool = False):
        """Initialize the SculptOptimizer with custom restraints.
        
        Args:
            custom_bonds: List of custom bond restraints to apply during optimization.
            custom_angles: List of custom angle restraints to apply during optimization.
            custom_torsions: List of custom torsion restraints to apply during optimization.
            full_sim: Not implemented. If True, run a full simulation instead of just minimization (1 ns). This can be useful to "shakeout" the structure to better satisfy difficult restraints.
        """
        self.custom_bonds = custom_bonds
        self.custom_angles = custom_angles
        self.custom_torsions = custom_torsions
        self.full_sim = full_sim  # By default, just do minimization

    def optimize(self, design: Design, sdf_files: List[str] = [], unique_naming: bool = False) -> Design:
        """Use EasyMD to optimize the structure file using custom restraints.
        
        Args:
            design: The Design structure file to optimize. Must have a 'structure' attribute.
            sdf_files: A list of SDF files to use for ligand parameterization.
            unique_naming: If True, generate a unique name for the output Design.
        Returns:
            A new Design object with the optimized structure.
        """
        # Create a temp directory:
        with tempfile.TemporaryDirectory() as temp_dir, \
             design.temp_structure_file() as temp_structure_file:

            ribbon.EasyMD(
                input_file=temp_structure_file,
                output_prefix=str(Path(temp_dir) / "optimized"),
                ligand_files=sdf_files,
                duration=1,
                custom_bonds=self.custom_bonds,
                custom_torsions=self.custom_torsions,
                custom_angles=self.custom_angles,
                minimize_only= True, # Eventually, add the full_sim flag. This will require me to pull coords in from a DCD.
            ).run()

            # What's the file that has our optimized PDB?
            output_pdb_file = Path(temp_dir) / "optimized_EM.pdb"

            # Create the output Design object:
            if unique_naming:
                output_name = generate_id()
            else:
                output_name = design.name + "_opt"

            # Make the new Design object:
            optimized_design = design.copy(set_parent=True)
            optimized_design.name = output_name
            optimized_design.load_structure(output_pdb_file)

            # Standardize the structure to change auth_seq_id to label_seq_id:
            standard_molecules = [StandardMolecule(structure_file = sdf) for sdf in sdf_files]
            optimized_design.structure.standardize(standard_molecules=standard_molecules, use_hydrogens=True, renumber=True)

            # Add to the history:
            optimized_design.history.append("Optimized using EasyMD. \n"
                                            "Custom bonds: {}, \n"
                                            "Custom angles: {}, \n"
                                            "Custom torsions: {}".format(
                self.custom_bonds, self.custom_angles, self.custom_torsions
            ))

        return optimized_design

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



    
