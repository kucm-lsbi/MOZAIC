from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski
from collections import namedtuple

# ESOL
class ESOLCalculator:
    def __init__(self):
        self.aromatic_query = Chem.MolFromSmarts("a")
        self.Descriptor = namedtuple("Descriptor", "mw logp rotors ap")

    def calc_ap(self, mol):
        matches = mol.GetSubstructMatches(self.aromatic_query)
        return len(matches) / mol.GetNumAtoms()

    def calc_esol_descriptors(self, mol):
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        rotors = Lipinski.NumRotatableBonds(mol)
        ap = self.calc_ap(mol)
        return self.Descriptor(mw=mw, logp=logp, rotors=rotors, ap=ap)

    def calc_esol(self, mol):
        intercept = 0.26121066137801696
        coef = {
            'mw': -0.0066138847738667125,
            'logp': -0.7416739523408995,
            'rotors': 0.003451545565957996,
            'ap': -0.42624840441316975
        }
        desc = self.calc_esol_descriptors(mol)
        esol = intercept + coef["logp"] * desc.logp + coef["mw"] * desc.mw + \
               coef["rotors"] * desc.rotors + coef["ap"] * desc.ap
        return esol
        
esol_calculator = ESOLCalculator()