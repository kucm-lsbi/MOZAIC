import logging
from vina import Vina
from .prep_docking import prepare_ligand

def run_vina(smiles, receptor_path, binding_center, docking_cfg):

    # Docking configs
    raw_box_size = docking_cfg['box_size']
    exhaustiveness = docking_cfg['exhaustiveness']
    n_poses = docking_cfg['n_poses']

    # Ligand
    ligand_pdbqt = prepare_ligand(smiles)
    
    if ligand_pdbqt is None:
        return None, None
        
    final_center = [float(c) for c in binding_center]
    final_box_size = [float(b) for b in raw_box_size]

    # Vina
    try:
        v = Vina(sf_name='vina', verbosity=0, cpu=1)
        v.set_receptor(receptor_path)
        v.set_ligand_from_string(ligand_pdbqt)
        
        v.compute_vina_maps(center=final_center, box_size=final_box_size)
        v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)

        energies = v.energies(n_poses=1)
        if energies is None or len(energies) == 0:
            logging.warning(f"No poses found for {smiles}")
            return None, None
            
        best_affinity = float(energies[0][0])
        
        pose_str = v.poses(n_poses=1, coordinates_only=False)
        
        return best_affinity, pose_str

    except Exception as e:
        logging.error(f"Vina docking failed for {smiles}: {e}")
        return None, None