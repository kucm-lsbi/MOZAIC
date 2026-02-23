import random
import copy
import logging

from rdkit import Chem

from ..utils.mol_check import validate_mol
from ..bank.mol_info import add_unique_ids
from ..utils.library import FUNCTIONAL_GROUPS, RXN_RULES, SMILES_FRAGMENTS
from .update_smi import update_smi

def get_compatible_fgs(initial_mol, target_uids):
    compatible_fgs = []
    
    target_set = set()
    raw_type = "Unknown"
    
    for item in target_uids:
        try:
            if isinstance(item, (list, tuple)):
                target_set.add(int(item[1]))
                raw_type = "Tuple/List"
            elif isinstance(item, (int, str)):
                target_set.add(int(item))
                raw_type = "Int/Str"
        except:
            continue
            
    for fg_name, smarts in FUNCTIONAL_GROUPS.items():
        pattern = Chem.MolFromSmarts(smarts)
        if not pattern: continue
        
        matches = initial_mol.GetSubstructMatches(pattern)
        
        for match in matches:
            match_uids = []
            for idx in match:
                atom = initial_mol.GetAtomWithIdx(idx)
                if atom.HasProp('unique_id'):
                    match_uids.append(int(atom.GetProp('unique_id')))
            
            if target_set.issubset(set(match_uids)):
                compatible_fgs.append(fg_name)
                break
       
    return compatible_fgs

def perform_mutation(item, initial_smiles, initial_mol, existing_smiles):

    new_item = copy.deepcopy(item)
    
    if not new_item['rxn_history']:
        return None

    entry = random.choice(new_item['rxn_history'])
    used_atoms = entry.get('used_atoms', [])

    # Skip
    if random.random() < 0.1:
        entry['selected_fragment'] = 'skip'
        
    else:
        # FG re-selection
        compatible_fgs = get_compatible_fgs(initial_mol, used_atoms)
            
        selected_fg = random.choice(compatible_fgs)
        entry['functional_group'] = selected_fg

        # Position re-selection
        compatible_rules = [
            r_name for r_name, r_info in RXN_RULES.items() 
            if selected_fg in r_info['functional_groups']
        ]
            
        selected_rule_name = random.choice(compatible_rules)
        rule_info = RXN_RULES[selected_rule_name]
        
        # Fragment re-selection
        if rule_info['num_reactants'] == 1:
            entry['selected_fragment'] = 'solo'
            entry['rxn_position'] = f"{selected_rule_name}_0"
            
        else:
            try:
                fg_idx = rule_info['functional_groups'].index(selected_fg)
            except ValueError:
                return None
            
            # My position
            current_pos_idx = str(fg_idx) # "0" or "1"
            entry['rxn_position'] = f"{selected_rule_name}_{current_pos_idx}"
            
            # Fragment position
            target_pos_idx = "1" if current_pos_idx == "0" else "0"
            compatible_frag_key = f"{selected_rule_name}_{target_pos_idx}"
            
            if compatible_frag_key not in SMILES_FRAGMENTS:
                return None
            
            frag_list = SMILES_FRAGMENTS[compatible_frag_key]
            if not frag_list:
                return None
                
            entry['selected_fragment'] = random.choice(frag_list)

    # Update SMILES
    prod = update_smi(initial_smiles, new_item)

    if prod is None:
        return None

    if prod == initial_smiles:
        return None

    if prod in existing_smiles:
        return None

    if validate_mol(prod) is None:
        return None

    new_item['product_smiles'] = prod
    new_item['scores'] = {}
    
    return new_item

def run_mutation(crossover_population, backup_crossover, seed, initial_smiles):

    target_size = len(seed) * 3
    backup_size = len(seed) * 3
    
    # succeeding crossover
    new_gen = list(crossover_population)
    backup = list(backup_crossover)
    
    existing_smiles = set()
    for item in new_gen + backup:
        if 'product_smiles' in item:
            existing_smiles.add(item['product_smiles'])
            
    raw_mol = Chem.MolFromSmiles(initial_smiles)

    initial_mol = add_unique_ids(raw_mol)
    
    # Filling
    def fill_list(target_list, limit, source_seeds):
        attempts = 0
        max_attempts = limit * 300
        
        while len(target_list) < limit and attempts < max_attempts:
            attempts += 1
            
            parent = random.choice(source_seeds)

            # initial_mol (with IDs)
            mutant = perform_mutation(parent, initial_smiles, initial_mol, existing_smiles)
            
            if mutant:
                target_list.append(mutant)
                existing_smiles.add(mutant['product_smiles'])


    # New_gen
    fill_list(new_gen, target_size, seed)
    
    # Backup
    fill_list(backup, backup_size, seed)
    
    logging.info(f"Mutation Complete. Final Sizes -> NewGen: {len(new_gen)}, Backup: {len(backup)}\n")

    return new_gen, backup