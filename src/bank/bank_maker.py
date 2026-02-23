import copy
import logging
import multiprocessing as mp

from rdkit import rdBase
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.SimDivFilters import rdSimDivPickers
from rdkit.Chem import rdFingerprintGenerator

from .growing import run_growing
from .mol_info import add_unique_ids
from ..utils.mol_check import validate_mol
from ..scoring.scorer import get_scores

def _worker_generate_chunk(args):

    initial_smiles, active_fgs_atoms, quota, worker_id = args
    
    rdBase.DisableLog('rdApp.error')
    rdBase.DisableLog('rdApp.warning')
    
    # Local Storage
    local_pool = []
    local_smiles = set()
    
    # Initial History Template
    initial_rxn_info = [{
        "step": 0, 
        "functional_group": fg, 
        "used_atoms": [], 
        "rxn_position": '', 
        "selected_fragment": ''
    } for fg in active_fgs_atoms.keys()]

    # Monitoring Variables
    total_attempts = 0
    failure_count = 0

    monitor_threshold = quota 

    # Generation Loop
    while len(local_pool) < quota:
        total_attempts += 1
        
        # growing
        product_smiles, product_rxn_info, inactive_atoms = run_growing(
            initial_smiles, initial_rxn_info, active_fgs_atoms
        )
        
        is_valid = False
        
        # product validation 
        if product_smiles is not None:
            if product_smiles not in local_smiles:
                if validate_mol(product_smiles) is not None:
                    all_initial = set(atom[1] for fg in active_fgs_atoms.values() for atoms in fg for atom in atoms)
                    used = set(atom[1] for atom in inactive_atoms)
                    
                    if all_initial.issubset(used):
                        is_valid = True
                        local_pool.append({
                            "product_smiles": product_smiles,
                            "rxn_history": product_rxn_info,
                            "scores": {}
                        })
                        local_smiles.add(product_smiles)
        
        # Failure Monitoring
        if not is_valid:
            failure_count += 1

        if total_attempts == monitor_threshold:
            fail_rate = (failure_count / monitor_threshold) * 100
            
            logging.info(f" >> [Worker {worker_id}] First {monitor_threshold} attempts: {fail_rate:.1f}% Invalid/Duplicate.\n")
            
            if fail_rate >= 95.0:
                error_msg = (f"  ❗️  Critical Failure Rate ({fail_rate:.1f}%).\n")
                logging.error(error_msg)
                raise RuntimeError(error_msg)

                
    return local_pool
    
def generate_diverse_bank(initial_smiles, active_fgs_atoms, n_bank, n_gen0, n_jobs):
    
    logging.info(f"Generating Gen0 pool ({n_gen0}) using {n_jobs} CPUs...\n")

    # Quota
    quotas = [n_gen0 // n_jobs] * n_jobs
    quotas[-1] += n_gen0 % n_jobs
    
    tasks = [
        (initial_smiles, active_fgs_atoms, q, i) 
        for i, q in enumerate(quotas)
    ]

    raw_results = []
    with mp.Pool(n_jobs) as pool:
        try:
            results_list = pool.map(_worker_generate_chunk, tasks)
            for res in results_list:
                raw_results.extend(res)
        except RuntimeError as e:
            raise e

    unique_pool = {}
    for item in raw_results:
        smi = item['product_smiles']
        if smi not in unique_pool:
            mol = Chem.MolFromSmiles(smi) 
            if mol:
                item['mol_obj'] = mol
                unique_pool[smi] = item
    
    generation_0 = list(unique_pool.values())
    logging.info(f"  > Generated {len(generation_0)} unique candidates (Target: {n_gen0}).\n")

    if len(generation_0) < n_bank:
        raise ValueError(f"Generated pool ({len(generation_0)}) is smaller than Bank size ({n_bank}). Check constraints.")

    # MaxMinPicker
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=4, fpSize=2048)
    fps = [mfpgen.GetFingerprint(item['mol_obj']) for item in generation_0]
    
    picker = rdSimDivPickers.MaxMinPicker()
    
    n_pick = min(len(generation_0), n_bank * 2)
    
    try:
        selected_indices = picker.LazyBitVectorPick(fps, len(fps), n_pick, seed=42)
    except AttributeError:
        def dist_func(i, j): return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
        selected_indices = picker.LazyPick(dist_func, len(fps), n_pick, seed=42)
        
    selected_indices = list(selected_indices)

    final_pool = []
    for idx in selected_indices:
        item = generation_0[idx]
        del item['mol_obj']
        final_pool.append(item)

    bank = final_pool[:n_bank]
    backup_bank = final_pool[n_bank:]
    
    print(f"Bank selection complete. (Bank: {len(bank)}, Backup: {len(backup_bank)})\n")
    
    return bank, backup_bank

def make_bank(initial_smiles, active_fgs_atoms, receptor_path, binding_center, docking_cfg, csa_cfg):

    n_bank   = csa_cfg['n_bank']
    n_gen0   = csa_cfg['n_gen0']
    n_jobs   = docking_cfg['n_cpu']
    
    # Make Bank
    bank, backup_bank = generate_diverse_bank(initial_smiles, active_fgs_atoms, n_bank, n_gen0, n_jobs)

    # Scoring
    scored_bank = get_scores(
        bank, 
        backup_bank, 
        receptor_path, 
        binding_center, 
        docking_cfg
    )

    return scored_bank