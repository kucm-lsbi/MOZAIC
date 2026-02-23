import random
import copy
import logging

from rdkit import Chem

from ..utils.mol_check import validate_mol
from .update_smi import update_smi

def get_matching_candidates(seeds, bank):

    candidates = []

    for s_idx, seed in enumerate(seeds):
        for b_idx, bank_item in enumerate(bank):
            
            for s_step_idx, s_step in enumerate(seed['rxn_history']):
                for b_step_idx, b_step in enumerate(bank_item['rxn_history']):
                    
                    # Same used atoms
                    s_atoms = sorted(s_step['used_atoms'])
                    b_atoms = sorted(b_step['used_atoms'])
                    
                    if s_atoms == b_atoms:
                        candidates.append((s_idx, b_idx, s_step_idx, b_step_idx))
                        
    return candidates

def execute_swap(seed_item, bank_item, s_step_idx, b_step_idx, initial_smiles, existing_smiles):

    child_seed = copy.deepcopy(seed_item)
    child_bank = copy.deepcopy(bank_item)
    
    e1 = child_seed['rxn_history'][s_step_idx]
    e2 = child_bank['rxn_history'][b_step_idx]
    
    # Swap
    e1['functional_group'],   e2['functional_group']   = e2['functional_group'],   e1['functional_group']
    e1['rxn_position'],       e2['rxn_position']       = e2['rxn_position'],       e1['rxn_position']
    e1['selected_fragment'],  e2['selected_fragment']  = e2['selected_fragment'],  e1['selected_fragment']

    prod1 = update_smi(initial_smiles, child_seed)
    prod2 = update_smi(initial_smiles, child_bank)
    
    result_seed = None
    result_bank = None
    
    # Validate Child 1
    if prod1 and prod1 != initial_smiles and prod1 not in existing_smiles:
        if validate_mol(prod1):
            child_seed['product_smiles'] = prod1
            child_seed['scores'] = {}
            result_seed = child_seed

    # Validate Child 2
    if prod2 and prod2 != initial_smiles and prod2 not in existing_smiles:
        if validate_mol(prod2):
            child_bank['product_smiles'] = prod2
            child_bank['scores'] = {}
            result_bank = child_bank
            
    return result_seed, result_bank
    
def calculate_quotas(candidates_init, candidates_curr, target_total):

    target_half = target_total // 2
    
    len_init = len(candidates_init)
    len_curr = len(candidates_curr)

    # Basic
    count_init = min(len_init, target_half)
    count_curr = min(len_curr, target_half)
    
    # Cross-filling
    if count_init < target_half:
        deficit = target_half - count_init
        count_curr = min(len_curr, count_curr + deficit)
        
    elif count_curr < target_half:
        deficit = target_half - count_curr
        count_init = min(len_init, count_init + deficit)
        
    return count_init, count_curr
    
def process_batch(candidate_list, seed, source_bank, dest_list, initial_smiles, existing_smiles, max_needed):

    added_count = 0
    
    for s_idx, b_idx, s_step, b_step in candidate_list:
        if len(dest_list) >= max_needed:
            break
            
        seed_obj = seed[s_idx]
        bank_obj = source_bank[b_idx]
        
        child1, child2 = execute_swap(
            seed_obj, bank_obj, s_step, b_step, 
            initial_smiles, existing_smiles
        )
        
        if child1:
            existing_smiles.add(child1['product_smiles'])
            dest_list.append(child1)
            added_count += 1
            
        if child2 and len(dest_list) < max_needed:
            existing_smiles.add(child2['product_smiles'])
            dest_list.append(child2)
            added_count += 1
            
    return added_count
    
def run_crossover(seed, initial_bank, current_bank, initial_smiles):
    
    new_gen = []
    backup = []
    existing_smiles = set()
    
    for s in seed: existing_smiles.add(s['product_smiles'])
    for b in initial_bank: existing_smiles.add(b['product_smiles'])
    for b in current_bank: existing_smiles.add(b['product_smiles'])
        
    # Matching candidates
    initial_candidates = get_matching_candidates(seed, initial_bank)
    current_candidates = get_matching_candidates(seed, current_bank)
    
    random.shuffle(initial_candidates)
    random.shuffle(current_candidates)
    
    target_gen_size = len(seed) * 2
    target_backup_size = len(seed) * 2

    # Generation
    n_init, n_curr = calculate_quotas(initial_candidates, current_candidates, target_gen_size)
    
    batch_init = initial_candidates[:n_init]
    batch_curr = current_candidates[:n_curr]
    
    # Execute swap
    process_batch(batch_init, seed, initial_bank, new_gen, initial_smiles, existing_smiles, target_gen_size)
    process_batch(batch_curr, seed, current_bank, new_gen, initial_smiles, existing_smiles, target_gen_size)
    
    # Del used
    remaining_init = initial_candidates[n_init:]
    remaining_curr = current_candidates[n_curr:]
    
    # Backup
    n_init_bk, n_curr_bk = calculate_quotas(remaining_init, remaining_curr, target_backup_size)
    
    batch_init_bk = remaining_init[:n_init_bk]
    batch_curr_bk = remaining_curr[:n_curr_bk]
    
    process_batch(batch_init_bk, seed, initial_bank, backup, initial_smiles, existing_smiles, target_backup_size)
    process_batch(batch_curr_bk, seed, current_bank, backup, initial_smiles, existing_smiles, target_backup_size)
    
    logging.info(f"Crossover Result -> NewGen: {len(new_gen)}, Backup: {len(backup)}")

    return new_gen, backup
