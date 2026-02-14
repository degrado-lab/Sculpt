from shim.fix import fix_atom_names_in_residue 
from pathlib import Path

input_dir = Path('./data/5RGA/')

input_file = input_dir / '5RGA_KEMP1_input.cif'
sdf_file = input_dir / 'KEMP1_h.sdf'

output_file = input_dir / '5RGA_KEMP1_fixed.cif'

fix_atom_names_in_residue(input_file, str(output_file), 'LIG2', sdf_file, use_hydrogens=False)
