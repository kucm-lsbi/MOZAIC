from rdkit import Chem
from rdkit.Chem import QED

def qed(smiles: str):

    mol = Chem.MolFromSmiles(smiles)

    return QED.qed(mol)
