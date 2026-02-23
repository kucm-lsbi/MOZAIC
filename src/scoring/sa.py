from rdkit import Chem
from rdkit.Contrib.SA_Score import sascorer

def sa(smiles: str):
    
    mol = Chem.MolFromSmiles(smiles)
    
    return sascorer.calculateScore(mol)
    