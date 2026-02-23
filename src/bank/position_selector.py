import random
import logging
from rdkit import Chem
from rdkit.Chem import rdFMCS

from ..utils.library import SMILES_FRAGMENTS, FUNCTIONAL_GROUPS, RXN_RULES
from .mol_info import add_unique_ids

def mapping_mcs(current_smiles, initial_smiles):
    
    current_mol = Chem.MolFromSmiles(current_smiles)
    initial_mol = Chem.MolFromSmiles(initial_smiles)
    
    if current_mol is None or initial_mol is None:
        return None, None
        
    # Initial molecule atom index mapping
    initial_mol = add_unique_ids(initial_mol)

    # MCS
    mcs_result = rdFMCS.FindMCS(
        [current_mol, initial_mol],
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        ringMatchesRingOnly=True, 
        completeRingsOnly=True
    )
    
    if not mcs_result or not mcs_result.smartsString:
        return None, None
    
    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
    if mcs_mol is None or mcs_mol.GetNumAtoms() == 0:
        return None, None
        
    # Initial -> MCS ID mapping
    match_initial = initial_mol.GetSubstructMatch(mcs_mol)
    if not match_initial:
        return None, None
        
    for mcs_idx, init_idx in enumerate(match_initial):
        init_atom = initial_mol.GetAtomWithIdx(init_idx)
        if init_atom.HasProp('unique_id'):
            mcs_mol.GetAtomWithIdx(mcs_idx).SetProp('unique_id', init_atom.GetProp('unique_id'))

    return Chem.MolToSmiles(mcs_mol), mcs_mol

def select_position(current_smiles, initial_smiles, available_fg, used_positions):

    mol = Chem.MolFromSmiles(current_smiles)

    # Get mapped mcs
    mcs_smiles, mcs_mol = mapping_mcs(current_smiles, initial_smiles)
    if mcs_mol:
        match = mol.GetSubstructMatch(mcs_mol)
        for mcs_idx, mol_idx in enumerate(match):
            mcs_atom = mcs_mol.GetAtomWithIdx(mcs_idx)
            if mcs_atom.HasProp('unique_id'):
                mol.GetAtomWithIdx(mol_idx).SetProp('unique_id', mcs_atom.GetProp('unique_id'))
        
    # Random functinoal group selection
    while available_fg:
        selected_fg = random.choice(list(available_fg.keys()))
        atoms_list = available_fg[selected_fg]

        functional_groups_in_mol = []
        smarts = FUNCTIONAL_GROUPS[selected_fg]
        pattern = Chem.MolFromSmarts(smarts)

        if mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            for match in matches:
                # matched atoms ids
                atoms_with_ids = []
                for idx in match:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.HasProp('unique_id'):
                        atoms_with_ids.append((atom.GetSymbol(), atom.GetProp('unique_id')))
                
                # check matched atoms
                if len(atoms_with_ids) == len(match):
                    if atoms_with_ids in atoms_list:
                         functional_groups_in_mol.append({
                                'functional_group': selected_fg, 
                                'atoms_with_ids': atoms_with_ids
                            })

        if not functional_groups_in_mol:
            # delete used fg
            available_fg.pop(selected_fg)
            continue

        fg_info = random.choice(functional_groups_in_mol)
        fg_key = fg_info['functional_group']
        atoms_with_ids = fg_info['atoms_with_ids']

        # Random position selection
        reactions_and_positions = []
        for rxn_name, rxn_info in RXN_RULES.items():
            if fg_key in rxn_info['functional_groups']:
                
                position = rxn_info['functional_groups'].index(fg_key)
                rxn_location = f"{rxn_name}_{position}"

                # used check
                if rxn_location not in used_positions.get(fg_key, []):
                    reactions_and_positions.append(rxn_location)
                    
        if not reactions_and_positions:
            # delete used fg
            available_fg.pop(selected_fg)
            continue
            
        rxn_position = random.choice(reactions_and_positions)
        
        rxn_info = [{
            "functional_group": fg_key, 
            "rxn_position": rxn_position, 
            "selected_fragment": '', 
            "step": 0, 
            "used_atoms": atoms_with_ids
        }]

        return rxn_position, rxn_info, fg_key, atoms_with_ids, available_fg

    # all used
    return None, None, None, None, available_fg

def update_used_positions(used_positions, inactive_atoms, available_fg, fg_key, rxn_position, used_atoms):

    # Used position
    if fg_key not in used_positions:
        used_positions[fg_key] = []
    used_positions[fg_key].append(rxn_position)

    # Used atoms
    for atom in used_atoms:
        if not any(ia[1] == atom[1] for ia in inactive_atoms):
            inactive_atoms.append(atom)

    # Update available fg
    used_uids = set(atom[1] for atom in inactive_atoms)
    
    keys_to_remove = []
    for fg_name, instances in available_fg.items():
        # survive atoms
        alive_instances = []
        for instance in instances:
            if not any(atom[1] in used_uids for atom in instance):
                alive_instances.append(instance)
        
        available_fg[fg_name] = alive_instances
        
        if not alive_instances:
            keys_to_remove.append(fg_name)
            
    for k in keys_to_remove:
        del available_fg[k]