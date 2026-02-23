import os

def parse_residues(site_list):

    # ['A:109,110,113', 'B:210,211'] --> {('A', 109), ('A', 110), ...}
    
    target_residues = set()

    for segment in site_list:
        try:
            chain_id, residues_str = segment.split(":", 1)
            chain_id = chain_id.strip()
            
            for res in residues_str.split(","):
                res_num = int(res.strip())
                target_residues.add((chain_id, res_num))
                
        except ValueError:
            raise ValueError(f"Invalid site format: '{segment}'. Expected 'Chain:Residues' (e.g., A:100,101)")

    return target_residues

def calculate_geometric_center(receptor_path, target_residues):

    atom_coords = []
    
    # PDBQT format
    # Chain ID: 21, ResSeq: 22-26, X: 30-38, Y: 38-46, Z: 46-54
    
    with open(receptor_path, 'r') as file:
        for line in file:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    chain_id = line[21].strip()
                    residue_index = int(line[22:26].strip())

                    if (chain_id, residue_index) in target_residues:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        atom_coords.append((x, y, z))
                        
                except (ValueError, IndexError):
                    continue

    if not atom_coords:
        return None

    # x_avg, y_avg, z_avg
    center = [sum(coord) / len(atom_coords) for coord in zip(*atom_coords)]
    
    return [round(c, 3) for c in center]

def prepare_center(receptor_path: str, site_list: list):

    parsed_set = parse_residues(site_list)
    
    center = calculate_geometric_center(receptor_path, parsed_set)
    
    if center is None:
        raise ValueError(
            f"Failed to compute docking center.\n"
            f"No atoms found in '{os.path.basename(receptor_path)}' matching sites: {site_list}\n"
            f"Check if Chain ID and Residue Numbers exist in the PDBQT file."
        )
        
    return center