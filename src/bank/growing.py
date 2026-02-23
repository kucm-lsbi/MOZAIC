import random

from rdkit import rdBase

rdBase.DisableLog('rdApp.error') 
rdBase.DisableLog('rdApp.warning')

from rdkit import Chem
from rdkit.Chem import AllChem

from ..utils.library import SMILES_FRAGMENTS, RXN_RULES
from .position_selector import mapping_mcs, select_position, update_used_positions

def sanitize_product(mol):

    try:
        Chem.SanitizeMol(mol)

        smi = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        return smi
    except Exception:
        return None
        
def run_growing(initial_smiles, initial_rxn_info, active_fgs_atoms):
    
    current_smiles = initial_smiles
    
    product_rxn_info = [info.copy() for info in initial_rxn_info]

    used_positions = {}
    inactive_atoms = []
    available_fg = active_fgs_atoms.copy()
    step = 1

    # Terminal check setup
    all_initial_atoms = set()
    for fg_atoms in active_fgs_atoms.values():
        for atoms in fg_atoms:
            for atom in atoms:
                all_initial_atoms.add(atom[1]) # Unique ID only
                
    while True:
        if not available_fg:
            break
            
        rxn_position, rxn_info, fg_key, used_atoms, available_fg = select_position(
            current_smiles, initial_smiles, available_fg, used_positions
        )
        
        if rxn_position is None:
            break

        # 10% skip
        if random.random() < 0.1:
            update_used_positions(used_positions, inactive_atoms, available_fg, fg_key, rxn_position, used_atoms)
            
            product_rxn_info.append({
                "step": step,
                "functional_group": fg_key,
                "used_atoms": used_atoms,
                "rxn_position": rxn_position,
                "selected_fragment": 'skip',
                "reacted_compound" : current_smiles
            })
            step += 1
            continue
            
        rxn_name, pos_idx = rxn_position.rsplit('_', 1)
        
        rule = RXN_RULES[rxn_name]
        num_reactants = rule['num_reactants']
        
        fragment_mol = None
        selected_frag_label = ''

        # Need fragment or not
        if num_reactants == 2:
            target_pos_idx = "1" if pos_idx == "0" else "0"
            compatible_frag_key = f"{rxn_name}_{target_pos_idx}"

            fragment_smiles = random.choice(SMILES_FRAGMENTS[compatible_frag_key])
            fragment_mol = Chem.MolFromSmiles(fragment_smiles)

            if fragment_mol is None:
                continue

            selected_frag_label = fragment_smiles

        elif num_reactants == 1:
            fragment_mol = None
            selected_frag_label = 'solo'

        current_step_info = rxn_info[0]
        current_step_info["selected_fragment"] = selected_frag_label
        current_step_info["step"] = step
        current_step_info["used_atoms"] = used_atoms
        
        # Reaction
        reaction_smarts = rule['reaction_string']
        rdkit_rxn = AllChem.ReactionFromSmarts(reaction_smarts)
        
        current_mol = Chem.MolFromSmiles(current_smiles)
        if current_mol is None: 
            continue

        # Initial atom index mapping
        mcs_smi, mcs_mol = mapping_mcs(current_smiles, initial_smiles)
        if mcs_mol:
            match = current_mol.GetSubstructMatch(mcs_mol)
            for i, idx in enumerate(match):
                atom_in_mcs = mcs_mol.GetAtomWithIdx(i)
                if atom_in_mcs.HasProp('unique_id'):
                    current_mol.GetAtomWithIdx(idx).SetProp('unique_id', atom_in_mcs.GetProp('unique_id'))

        products_list = []
        try:
            if num_reactants == 2:
                if pos_idx == "0":
                    products_list = rdkit_rxn.RunReactants((current_mol, fragment_mol))
                else:
                    products_list = rdkit_rxn.RunReactants((fragment_mol, current_mol))
            elif num_reactants == 1:
                products_list = rdkit_rxn.RunReactants((current_mol,))
        except Exception:
            products_list = []

        # Post-processing
        success = False

        for product_candidates in products_list:
            if not product_candidates:
                continue

            product_mol_obj = product_candidates[0]
            product_smiles_candidate = sanitize_product(product_mol_obj)

            if product_smiles_candidate:
                current_smiles = product_smiles_candidate
                success = True
                
                # History update
                product_rxn_info = [entry for entry in product_rxn_info if entry["step"] != 0]
                
                current_step_info["reacted_compound"] = current_smiles
                product_rxn_info.append(current_step_info)

                update_used_positions(used_positions, inactive_atoms, available_fg, fg_key, rxn_position, used_atoms)
                step += 1

                break

        if not success:
            continue
                    
        used_atom_ids = set(ia[1] for ia in inactive_atoms)
        if all_initial_atoms.issubset(used_atom_ids):
            return current_smiles, product_rxn_info, inactive_atoms

    return None, None, inactive_atoms