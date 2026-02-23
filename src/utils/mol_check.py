import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize as std

def auto_fix_quaternary_neutral_N(m):
    # Quaternary Neutral N to +1
    fixed = False
    for a in m.GetAtoms():
        if a.GetSymbol() == 'N': 
            val_exp = a.GetExplicitValence()
            if a.GetFormalCharge() == 0 and val_exp >= 4:
                a.SetFormalCharge(1)
                a.UpdatePropertyCache(strict=False)
                fixed = True
    return fixed

def validate_mol(smiles: str):

    # 1. smi to mol
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return None
    
    try:
        m.UpdatePropertyCache(strict=False)
    except Exception:
        return None

    # 2. Boron filtering (vina)
    if any(atom.GetAtomicNum() == 5 for atom in m.GetAtoms()):
        return None

    # 3. Standardazation
    try:
        m = std.Cleanup(m)
        m = std.Normalize(m)
        m = std.Reionize(m)
    except Exception:
        return None

    # 4. Quaternary Neutral N to +1
    auto_fix_quaternary_neutral_N(m)

    # 5. Full Sanitize Check
    try:
        Chem.SanitizeMol(m)
    except Exception:
        return None

    # 6. Round-trip Check
    try:
        smi2 = Chem.MolToSmiles(m, isomericSmiles=True)
        m2 = Chem.MolFromSmiles(smi2)
        Chem.SanitizeMol(m2)
    except Exception:
        return None

    # 7. Identifying weird atom or bond
    for a in m2.GetAtoms():
        if a.GetSymbol() == 'N' and a.GetFormalCharge() == 0:
            if a.GetExplicitValence() >= 4:
                return None
    for b in m2.GetBonds():
        if b.GetBondType() == Chem.BondType.UNSPECIFIED:
            return None

    # 8. 3D Embedding Test
    try:
        mh = Chem.AddHs(m2)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        params.maxIterations = 1000
        
        if AllChem.EmbedMolecule(mh, params) != 0:
            return None
    except Exception:
        return None
        
    return m2