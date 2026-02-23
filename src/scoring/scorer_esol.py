import multiprocessing as mp
import logging
from tqdm import tqdm

from rdkit import Chem  

from .qed import qed
from .sa import sa
from .docking import run_vina
from .esol import esol_calculator 

def calculate_chem_props(entry):

    smi = entry["product_smiles"]

    # QED, SA
    q = qed(smi)
    s = sa(smi)
    s_norm = (10.0 - s) / 9.0

    # ✅ ESOL
    mol = Chem.MolFromSmiles(smi)
    esol = esol_calculator.calc_esol(mol) if mol is not None else None

    entry.setdefault("scores", {})
    entry["scores"].update({
        "qed":      q,
        "sa":       s,
        "sa_norm":  s_norm,
        "esol":     esol,   # ✅ 추가
    })

def docking_worker(args):
    
    idx, entry, receptor_path, binding_center, docking_cfg = args
    smiles = entry["product_smiles"]
    
    aff, pose = run_vina(smiles, receptor_path, binding_center, docking_cfg)
    
    return idx, aff, pose
    
def get_scores(population, backup, receptor_path, binding_center, docking_cfg):

    n_jobs = docking_cfg['n_cpu']
    
    # QED, SA, ESOL
    for entry in population + backup:
        calculate_chem_props(entry)

    # Docking
    tasks_pop = [
        (i, entry, receptor_path, binding_center, docking_cfg)
        for i, entry in enumerate(population)
    ]

    failed_indices = []
    
    print(f"[Scoring] Docking {len(population)} candidates...\n")
    
    with mp.Pool(n_jobs) as pool:
        for idx, aff, pose in tqdm(pool.imap_unordered(docking_worker, tasks_pop), total=len(tasks_pop)):
            if aff is not None:
                population[idx]['scores']['affinity'] = aff
                population[idx]['scores']['pose'] = pose
            else:
                failed_indices.append(idx)

    # Back docking
    n_failures = len(failed_indices)
    
    if n_failures > 0:
        logging.warning(f"{n_failures} failures detected. Running backups...")
        
        tasks_backup = [
            (i, entry, receptor_path, binding_center, docking_cfg)
            for i, entry in enumerate(backup)
        ]
        
        successful_backups = []
        
        with mp.Pool(n_jobs) as pool:
            iterator = pool.imap_unordered(docking_worker, tasks_backup)
            
            for idx, aff, pose in tqdm(iterator, total=len(tasks_backup), desc="Docking Backups"):
                if aff is not None:
                    successful_backups.append((idx, aff, pose))

                    if len(successful_backups) >= n_failures:
                        break

        # Replace
        logging.info(f"Replacing {n_failures} failures with successful backups...")
        
        used_backup_count = 0
        for fail_idx in failed_indices:
            
            bk_idx, bk_aff, bk_pose = successful_backups.pop(0)
            replacement_entry = backup[bk_idx]
            
            population[fail_idx] = replacement_entry
            population[fail_idx]['scores']['affinity'] = bk_aff
            population[fail_idx]['scores']['pose'] = bk_pose
            
            used_backup_count += 1
            
        logging.info(f"Successfully replaced {used_backup_count}/{n_failures} failures.")

    return population