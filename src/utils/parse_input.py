import os
import yaml
import shutil
import logging

from openbabel import openbabel

obErrorLog = openbabel.OBMessageHandler()
obErrorLog.SetOutputLevel(0)

def load_config(config_path):

    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
        
    # Config to dict
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
        
    return config
    
def create_job_dir(output_root, job_name):

    if not os.path.exists(output_root):
        os.makedirs(output_root)

    base_path = os.path.join(output_root, job_name)
    target_path = base_path

    # Appends _1, _2 ... if exists
    counter = 1
    while os.path.exists(target_path):
        target_path = f"{base_path}_{counter}"
        counter += 1
    
    os.makedirs(target_path)
    
    return target_path

def validate_file(file_path, allowed_extensions=None):

    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Input file not found: {file_path}")
    
    # Check file extension
    if allowed_extensions:
        ext = os.path.splitext(file_path)[1].lower()
        allowed = [e.lower() for e in allowed_extensions]
        
        if ext not in allowed:
            raise ValueError(
                f"Invalid file extension: {file_path}\n"
                f"Allowed formats: {allowed_extensions}"
            )
            
    return os.path.abspath(file_path)

def load_smiles(input_smi):
    
    valid_path = validate_file(input_smi, allowed_extensions=['.smi'])
    
    with open(input_smi, 'r') as f:
        smiles = f.readline().strip()
        
        if not smiles:
            raise ValueError(f"SMILES file is empty: {input_smi}")
        return smiles

    # Raw
    return input_val.strip()

def pdb_to_pdbqt_rigid(input_path, output_path):

    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdbqt")
    
    obConversion.AddOption("r", openbabel.OBConversion.OUTOPTIONS) 

    mol = openbabel.OBMol()
    
    if not obConversion.ReadFile(mol, input_path):
        raise IOError(f"OpenBabel failed to read file: {input_path}")

    if not obConversion.WriteFile(mol, output_path):
        raise IOError(f"OpenBabel failed to write file: {output_path}")
        
def prepare_receptor(input_path, job_dir):

    valid_path = validate_file(input_path, allowed_extensions=['.pdb', '.pdbqt'])

    file_name = os.path.basename(valid_path)
    base_name, ext = os.path.splitext(file_name)
    ext = ext.lower()

    target_path = os.path.join(job_dir, base_name + ".pdbqt")
    
    try:
        if ext == '.pdbqt':
            logging.info(f"Copying receptor PDBQT: {valid_path} -> {target_path}\n")
            shutil.copy(valid_path, target_path)
            
        elif ext == '.pdb':
            logging.info(f"Converting receptor PDB to PDBQT: {valid_path} -> {target_path}\n")
            pdb_to_pdbqt_rigid(valid_path, target_path)
            
        return target_path

    except Exception as e:
        logging.error(f"Failed to prepare receptor: {e}")
        raise e