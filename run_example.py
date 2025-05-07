from sculpt import sculpt
from pathlib import Path

def main():
    input_dir = Path('./data/5RGA/')
    run_dir = Path('./5RGA_KEMp1/')

    # 🧪 Ligand information
    ligand_smiles = 'O=C1C=C(C)C(C=CC2=C3C=NO2)=C3O1'  # SMILES string for ligand (KEMp-1))
    ligand_sdf = Path('./data/5RGA/KEMP1_h.sdf')       # Path to ligand SDF file

    # 🧬 Protein redesign
    residues_to_change = "A125 A126 A127 A128 A130 A17 A170 A172 A19 A20 A205 A206 A207 A208 A209 A21 A22 A23 A234 A235 A236 A237 A238 A264 A265 A266 A267 A268 A401 A42 A43 A44 A45 A46 A47 A48 A49 A50 A502 A524 A551 A587 A611 A617 A633 A635 A669 A69 A692 A700 A730 A81 A82 A83 A84 A85 A86 A87 A90"
    # We'll allow several residues in the active site to change, but keeping residue 129 constant (the catalytic base).

    # 🔁 Pipeline parameters
    num_cycles = 10              # Number of Sculpt optimization cycles
    lmpnn_design_num  = 10        # Sequences to generate per structure (via LigandMPNN)
    top_structures_num = 4        # Keep this many top designs per cycle (x5 structures each)

    # 🖥️ Scheduler config
    scheduler = 'local'          # 'local', 'SGE', or 'SLURM'
    queue = '*'                  # Queue name (for job schedulers)
    node_name  = '*'              # Node name (for job schedulers)

    # 🎯 Optimization geometry constraints
    sculpt_distances = [    # Format: "CHAIN:RES:ATOM, CHAIN:RES:ATOM, force, target_dist"
        'A:129:OD1, B:1:H3, 20, 1.5'
    ]
    sculpt_torsions = [     # Format: "atom1, atom2, atom3, atom4, force, periodicity, target"
        'A:129:OD2, A:129:CG, A:129:OD1, B:1:H3, 20, 2, 0'
    ]
    sculpt_angles = [       # Format: "atom1, atom2, atom3, force, angle"
        'A:129:CG, A:129:OD1, B:1:H3, 50, 240',
        'A:129:OD1, B:1:H3, B:1:C10_, 50, 180'
    ]

    # 🧪 Scoring function
    # Currently, the scoring function is set to be the distance between the closest carboxylate oxygen in residue 129 and the ligand proton.
    # There's no user-interface for this currently - stay tuned!
    # If you want to change the scoring function, you can do so in the 'calculate_distance' function in sculpt.py.

    sculpt(
        input_dir=input_dir,
        run_dir=run_dir,
        num_cycles=num_cycles,
        ligand_smiles=ligand_smiles,
        ligand_SDF=ligand_sdf,
        residues_to_change=residues_to_change,
        lmpnn_design_num=lmpnn_design_num,
        top_structures_num=top_structures_num,
        scheduler=scheduler, # 'SGE' or 'SLURM' if you want to use a job scheduler.
        queue=queue,
        node_name=node_name,
        distances=sculpt_distances,
        torsions=sculpt_torsions,
        angles=sculpt_angles
    )

if __name__ == '__main__':
    main()
