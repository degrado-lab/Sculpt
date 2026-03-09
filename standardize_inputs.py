from sculpt import standardize

standardize(
    structure_file='run_data/data/5RGA/5RGA_KEMP1_input.cif',
    output_file='run_data/data/5RGA/5RGA_KEMP1_input_standardized_h.cif',
    sdf_files=['run_data/data/5RGA/KEMP1_h.sdf'],
    add_hydrogens = True
)

