import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs import TanimotoSimilarity

def distance(smiles1, smiles2):
    
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    # Radius 4, bit 2048 
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 4, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 4, nBits=2048)
    
    dist = 1 - TanimotoSimilarity(fp1, fp2)
    
    return dist

def calc_davg(generation):
    
    n = len(generation)
    
    total_distance = 0
    count = 0

    for i in range(n):
        for j in range(i + 1, n):
            d = distance(generation[i]['product_smiles'], generation[j]['product_smiles'])
            total_distance += d
            count += 1

    return total_distance / count if count > 0 else 0
    
def update_bank(current_bank, child, dcut):

    all_population = current_bank + child

    neg_affs = [-entry['scores']['affinity'] for entry in all_population]

    min_aff = min(neg_affs)
    max_aff = max(neg_affs)
    
    denom = max_aff - min_aff

    if denom < 1e-9:
        raise ValueError(
            f"Error: Max and Min affinities are identical ({min_aff}). "
        )
        
    for entry in all_population:
        aff = entry['scores']['affinity']
        qed = entry['scores']['qed']
        norm_sa = entry['scores']['sa_norm']

        neg_aff = -aff
        norm_aff = (neg_aff - min_aff) / denom
        entry['scores']['affinity_norm'] = norm_aff
        
        # Objective (QED + SA_norm + Vina_norm)
        score_metric = qed + norm_sa + norm_aff
        entry['scores']['objective'] = score_metric

    updated_gen = current_bank.copy()
    n_updated = 0

    # Update
    for cand in child:
        cand_score = cand['scores']['objective']

        dists = [distance(cand['product_smiles'], mol['product_smiles']) for mol in updated_gen]
        min_dist = min(dists)
        min_dist_idx = dists.index(min_dist)

        # Niche (1 vs 1)
        if min_dist < dcut:
            existing_score = updated_gen[min_dist_idx]['scores']['objective']
            if cand_score > existing_score:
                updated_gen[min_dist_idx] = cand
                n_updated += 1
                
        # Bank-wide (1 vs All)
        else:
            updated_gen.append(cand)
            
            all_scores = [mol['scores']['objective'] for mol in updated_gen]
            min_score_idx = np.argmin(all_scores)
            
            if min_score_idx != len(updated_gen) - 1:
                n_updated += 1
            
            updated_gen.pop(min_score_idx)

    return updated_gen, n_updated