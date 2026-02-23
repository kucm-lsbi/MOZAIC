from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from ..bank.position_selector import mapping_mcs, update_used_positions
from ..utils.library import RXN_RULES

def update_smi(initial_smiles, population):
    
    rdBase.DisableLog('rdApp.error')
    rdBase.DisableLog('rdApp.warning')
    
    current_smiles = initial_smiles
    
    used_positions = {}  
    inactive_atoms = []
    remaining_fgs = {}

    for step_info in population['rxn_history']:
        step_num          = step_info.get("step", 0)
        fg_key            = step_info.get("functional_group", "")
        used_atoms        = step_info.get("used_atoms", [])
        rxn_position      = step_info.get("rxn_position", "")
        selected_fragment = step_info.get("selected_fragment", None)

        if selected_fragment == 'skip' or not rxn_position:
            step_info['reacted_compound'] = current_smiles
            update_used_positions(used_positions, inactive_atoms, remaining_fgs, fg_key, rxn_position, used_atoms)
            continue

        # Reaction rules
        rxn_name = rxn_position.rsplit('_', 1)[0]
        reaction_smarts = RXN_RULES[rxn_name]['reaction_string']
        rxn = AllChem.ReactionFromSmarts(reaction_smarts)

        current_mol = Chem.MolFromSmiles(current_smiles)

        reactant_tuple = None
        core = None
        
        if selected_fragment == 'solo':
            reactant_tuple = (current_mol,)
            core = reactant_tuple[0]
            
        else:
            if not selected_fragment:
                return None
                
            fragment_mol = Chem.MolFromSmiles(selected_fragment)
            if fragment_mol is None: return None
                
            if rxn_position.endswith("_0"):
                reactant_tuple = (current_mol, fragment_mol)
                core = reactant_tuple[0]
            elif rxn_position.endswith("_1"):
                reactant_tuple = (fragment_mol, current_mol)
                core = reactant_tuple[1]

        if core is None:
            return None
            
        # Used atom mapping
        core_smiles = Chem.MolToSmiles(core)
        mcs_smiles_core, mcs_mol_core = mapping_mcs(core_smiles, initial_smiles)
        
        if mcs_smiles_core is not None and mcs_mol_core is not None:
            match = core.GetSubstructMatch(mcs_mol_core)
            for i, idx in enumerate(match):
                if mcs_mol_core.GetAtomWithIdx(i).HasProp('unique_id'):
                    core.GetAtomWithIdx(idx).SetProp('unique_id', mcs_mol_core.GetAtomWithIdx(i).GetProp('unique_id'))
        else:
            return None

        # Reaction
        try:
            products_list = rxn.RunReactants(reactant_tuple)
            if not products_list:
                    return None
    
            product_mol = products_list[0][0]
                
            Chem.SanitizeMol(product_mol)
            
            product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True, canonical=True)
            step_info['reacted_compound'] = product_smiles
            current_smiles = product_smiles
    
            update_used_positions(used_positions, inactive_atoms, remaining_fgs, fg_key, rxn_position, used_atoms)
            
        except Exception:
            return None
            
    return current_smiles
