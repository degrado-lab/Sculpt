from sculpt.design import Design
from sculpt.tasks.fold import SculptHydrogenAdder
from shim import StandardMolecule
from pathlib import Path

def standardize(
            structure_file,
            output_file,
            sdf_files: list = [],
            add_hydrogens = True,
            ):
    """
    Return a single input file, with standardized atom names and (optionally, but recommended) hydrogens added.
    This is useful for generating an initial input structure, 
    so you can reference the correct atom names.

    Args:
        structure_file: Path to the input structure file.
        sdf_files: List of paths to SDF files for standardization.
        add_hydrogens: Whether to add hydrogens to the structure.
    """

    # Grab the current input structures:
    file = Path(structure_file)
    if file.suffix in ['.pdb', '.cif', '.mmcif']:
        input_design = Design(name=file.stem)
        input_design.load_structure(file)
    else:
        raise ValueError(f'Unsupported file format: {file.suffix}')
    
    # Create standardized input structures:
    standard_molecules = [StandardMolecule(structure_file=sdf_file) for sdf_file in sdf_files]
    input_design.structure.standardize(standard_molecules=standard_molecules)#, use_hydrogens=False)
    

    # Add hydrogens if requested:
    if add_hydrogens:
        hydrogen_adder = SculptHydrogenAdder()
        input_design = hydrogen_adder.add_hydrogens(input_design, sdf_files)

    # Output the standardized structure:
    input_design.structure_file(output_file)
    