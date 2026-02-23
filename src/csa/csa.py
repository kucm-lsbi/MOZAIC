import os
import logging
from copy import deepcopy

from ..bank.mol_info import get_initial_mol
from ..bank import make_bank
from ..child import make_child
from ..utils.save_results import save_xlsx, save_pose

from .selection import make_seed, mark_unused, count_unused
from .update import calc_davg, update_bank

def run_csa(job_dir, initial_smiles, receptor_path, binding_center, docking_cfg, csa_cfg):
    
    active_fgs_atoms = get_initial_mol(job_dir, initial_smiles)
        
    logging.info("Generating initial bank...\n")   
    initial_bank = make_bank(
        initial_smiles,
        active_fgs_atoms,
        receptor_path,
        binding_center,
        docking_cfg,
        csa_cfg
    )
    
    print()
    save_xlsx(initial_bank, os.path.join(job_dir, "initial_bank.xlsx"))
    
    current_bank = deepcopy(initial_bank)
    current_bank = mark_unused(current_bank)
    
    n_bank = csa_cfg['n_bank']
    n_seed = int(n_bank * 0.5)
    
    # Dcut
    initial_davg = calc_davg(initial_bank)
    dcut = initial_davg / 2.0
    dcut_min = initial_davg / 5.0
    
    logging.info("-" * 25)
    logging.info(f"Initial Davg : {initial_davg:.3f}")
    logging.info(f"Start Dcut   : {dcut:.3f}")
    logging.info(f"Min Dcut     : {dcut_min:.3f}")
    logging.info(f"Seed Size    : {n_seed}")
    logging.info("-" * 25 + "\n")

    # Main loop
    n_rounds = csa_cfg['n_rounds']
    
    for round_num in range(1, n_rounds + 1):
        logging.info(f"[Round {round_num}/{n_rounds}] Started\n")
        
        while True:

            seed = make_seed(current_bank, n_seed)

            logging.info("Generating child ...\n")
            child = make_child(
                initial_bank,
                current_bank,
                seed,
                initial_smiles,
                receptor_path,
                binding_center,
                docking_cfg
            )

            child = mark_unused(child)
            
            # Bank Update
            current_bank, n_updated = update_bank(current_bank, child, dcut)
            
            # Unused Count
            unused_count = count_unused(current_bank)
            print()
            logging.info(f"  - Bank Updated: {n_updated} times")
            logging.info(f"  - Dcut: {dcut:.3f} | Unused: {unused_count}\n")

            # Dcut Annealing
            dcut = max(dcut_min, dcut * 0.98)

            # Round break
            if unused_count < n_seed:
                logging.info("  >> Unused threshold reached. Resetting used status.\n")
                current_bank = mark_unused(current_bank)
                break

    final_bank = current_bank
    save_xlsx(final_bank, os.path.join(job_dir, "final_bank.xlsx"))
    
    pose_dir = os.path.join(job_dir, "poses")
    if not os.path.exists(pose_dir):
        os.makedirs(pose_dir)
        
    save_pose(final_bank, pose_dir, receptor_path)
    
    return final_bank, initial_bank