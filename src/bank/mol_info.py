import sys
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from ..utils.library import FUNCTIONAL_GROUPS

def add_unique_ids(mol):

    for atom in mol.GetAtoms():
        atom.SetProp('unique_id', str(atom.GetIdx()))
        
    return mol

def find_functional_groups(mol):

    active_fgs = {}

    for fg_name, smarts in FUNCTIONAL_GROUPS.items():
        pattern = Chem.MolFromSmarts(smarts)
        
        if not mol.HasSubstructMatch(pattern):
            continue

        matches = mol.GetSubstructMatches(pattern)
        valid_matches = []

        for match in matches:
            # [Filter] Aromatic Nitrogen
            skip_match = False
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'N' and atom.GetIsAromatic():
                    skip_match = True
                    break
            
            if skip_match:
                continue
            
            # Extract Atom IDs
            atoms_with_ids = []
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                uid = atom.GetProp('unique_id') if atom.HasProp('unique_id') else str(idx)
                atoms_with_ids.append((atom.GetSymbol(), uid))
            
            valid_matches.append(atoms_with_ids)

        if valid_matches:
            active_fgs[fg_name] = valid_matches

    return active_fgs

def print_fg_analysis_table(active_fgs):

    rows = []
    for fg_name, matches in active_fgs.items():
        for atoms in matches:
            # atoms example: [('C', '2'), ('O', '3')]
            
            # [Main Atom Logic]
            main_atom = atoms[0]
            for sym, uid in atoms:
                if sym != 'C':
                    main_atom = (sym, uid)
                    break
            
            rows.append({
                'main_sym': main_atom[0],
                'main_uid': int(main_atom[1]),
                'fg_name': fg_name,
                'atoms': atoms
            })

    rows.sort(key=lambda x: (x['main_uid'], x['fg_name']))

    from collections import defaultdict
    grouped_rows = defaultdict(list)
    for row in rows:
        key = (row['main_sym'], row['main_uid'])
        grouped_rows[key].append(row)

    print("-"*100)
    print(f" {'No.':<4} | {'Main Atom (ID)':<15} | {'FG Name':<50} | {'Full Atoms'}")
    print("-" * 100)
    
    # Selection Map 
    selection_map = []
    idx = 1
    
    for (main_sym, main_uid), group in grouped_rows.items():
        
        group_data = []
        
        for i, row in enumerate(group):
            if i == 0:
                no_str = str(idx)
                main_str = f"{main_sym}:{main_uid}"
            else:
                no_str = ""
                main_str = ""
            
            full_str = ", ".join([f"{sym}:{uid}" for sym, uid in row['atoms']])
            display_name = (row['fg_name'][:47] + '..') if len(row['fg_name']) > 47 else row['fg_name']
            
            print(f" {no_str:<4} | {main_str:<15} | {display_name:<50} | [{full_str}]")
            
            group_data.append((row['fg_name'], row['atoms']))
        
        print("-" * 100 + "\n")
        
        selection_map.append(group_data)
        idx += 1
    
    return selection_map
    
def select_active_sites_interactive(active_fgs, img_path):

    # Table
    selection_map = print_fg_analysis_table(active_fgs)
    total_sites_count = len(selection_map)
    
    print(f"💡 Recommendation: Check the atom positions before selection '{img_path}'\n")
    
    print("Usage: Enter the 'No.' of sites to use (space separated, e.g. 1 2 3).\n")
    print("       Type 'all' to use everything, or 'q' to quit.\n")
    
    while True:
        user_input = input(" >> Select atoms to start MOZAIC: ").strip().lower()
        print()
        
        if not user_input:
            continue
            
        if user_input in ['q', 'quit', 'exit']:
            print(">> Operation cancelled.")
            sys.exit(0)
            
        selected_indices = []

        if user_input == 'all':
            selected_indices = list(range(1, total_sites_count + 1))
        else:
            try:
                selected_indices = [int(x) for x in user_input.split()]
                if any(i < 1 or i > total_sites_count for i in selected_indices):
                    print("Error: Invalid number selected. Check the table 'No.'")
                    continue
            except ValueError:
                print("Error: Please enter numbers only.")
                continue

        selected_count = len(selected_indices)
        
        if selected_count >= 3:
            print(f" ❗️ [WARNING] You selected {selected_count} active atoms.\n")
            print("  -- Growing 3 or more steps can raise the risk of generating invalid SMILES.\n")
            print("  -- MOZAIC will HALT if the generation failure rate exceeds 95%.\n")

        # 3. 딕셔너리 필터링 및 리턴
        selected_groups = [selection_map[i-1] for i in selected_indices]
        
        filtered_fgs = {}
        real_atom_count = 0
        
        for group in selected_groups:
            for sel_name, sel_atoms in group:
                if sel_name not in filtered_fgs:
                    filtered_fgs[sel_name] = []
                filtered_fgs[sel_name].append(sel_atoms)
            real_atom_count += 1
               
        print(f" >> Selected {real_atom_count} atoms.\n")
        return filtered_fgs
            
def visualize_initial_mol(mol, active_fgs, save_dir, filename="initial_mol.png"):

    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)
    full_path = save_path / filename

    rdDepictor.Compute2DCoords(mol)
    rdDepictor.NormalizeDepiction(mol)

    # Highlight active atoms
    highlight_atoms = set()
    for matches in active_fgs.values():
        for match in matches:
            for symbol, uid in match:
                highlight_atoms.add(int(uid))
    
    drawer = rdMolDraw2D.MolDraw2DCairo(800, 600)
    opts = drawer.drawOptions()
    opts.addAtomIndices = True          
    if hasattr(opts, "atomLabelFontSize"):
        opts.atomLabelFontSize = 20
    elif hasattr(opts, "minFontSize"):
        opts.minFontSize = 20     
    opts.setHighlightColour((1.0, 0.8, 0.0, 0.5))
    opts.highlightBondWidthMultiplier = 2

    drawer.DrawMolecule(mol, highlightAtoms=list(highlight_atoms))
    drawer.FinishDrawing()
    
    with open(full_path, 'wb') as f:
        f.write(drawer.GetDrawingText())
    
    print(f"[Visualization] Initial molecule saved to: {full_path}\n")
    return full_path

def get_initial_mol(job_dir, smiles):

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
        
    mol = add_unique_ids(mol)

    active_fgs = find_functional_groups(mol)
    if not active_fgs:
        raise ValueError("No reactive functional groups found in the initial smiles.")
        
    img_path = visualize_initial_mol(mol, active_fgs, job_dir)

    selected_fgs = select_active_sites_interactive(active_fgs, img_path)

    return selected_fgs