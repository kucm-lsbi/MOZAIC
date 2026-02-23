import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def smi_to_mol3d(smiles):

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)

    # Embedding
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    
    res = AllChem.EmbedMolecule(mol, params)
    
    if res == -1:
        params.useRandomCoords = True
        res = AllChem.EmbedMolecule(mol, params)
        if res == -1:
            logging.warning(f"3D Embedding Failed: {smiles}")
            return None

    # Optimization
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant='MMFF94s')
        else:
            AllChem.UFFOptimizeMolecule(mol)
            
    except Exception:
        pass
    
    return mol

def _extract_pdbqt_string(writer_out):

    if writer_out is None:
        return None

    if isinstance(writer_out, str):
        pdbqt = writer_out

    elif isinstance(writer_out, (bytes, bytearray)):
        pdbqt = writer_out.decode("utf-8", errors="replace")

    elif isinstance(writer_out, (tuple, list)):
        pdbqt = None

        first = writer_out[0] if len(writer_out) > 0 else None
        if isinstance(first, str):
            pdbqt = first
        elif isinstance(first, (list, tuple)) and all(isinstance(x, str) for x in first):
            pdbqt = "".join(first)

        if pdbqt is None:
            for item in writer_out:
                if isinstance(item, str):
                    pdbqt = item
                    break
                if isinstance(item, (list, tuple)) and all(isinstance(x, str) for x in item):
                    pdbqt = "".join(item)
                    break

        if pdbqt is None:
            return None

    else:
        return None

    if not any(tag in pdbqt for tag in ("ATOM", "HETATM", "ROOT", "TORSDOF")):
        logging.warning("Extracted PDBQT string does not look like a standard PDBQT block.")

    pdbqt = pdbqt.rstrip("\n") + "\n"
    return pdbqt


def mol3d_to_pdbqt(mol):
    try:
        prep = MoleculePreparation()
        mol_setups = prep.prepare(mol)

        if not mol_setups:
            logging.warning("Meeko conversion failed: empty mol_setups")
            return None

        writer_out = PDBQTWriterLegacy.write_string(mol_setups[0])

        pdbqt_string = _extract_pdbqt_string(writer_out)
        if pdbqt_string is None:
            logging.warning(f"Meeko conversion failed: unexpected writer output type={type(writer_out)}")
            return None

        return pdbqt_string

    except Exception as e:
        logging.warning(f"Meeko conversion failed: {e}")
        return None


# smi - 3d - pdbqt
def prepare_ligand(smiles):
    mol3d = smi_to_mol3d(smiles)
    if mol3d is None:
        return None

    return mol3d_to_pdbqt(mol3d)