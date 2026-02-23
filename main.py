import sys
import os
import argparse
import logging
import shutil

from datetime import datetime
from rdkit import Chem

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.warning')

from src.csa import run_csa 
from src.utils.logger import setup_logger
from src.utils.prepare_center import prepare_center
from src.utils.parse_input import (
    create_job_dir,
    prepare_receptor,
    load_config, 
    load_smiles
)

def parse_args():
    parser = argparse.ArgumentParser(
        description="MOZAIC: Molecule Optimization via SMARTS-based fragment growing & CSA",
        formatter_class=argparse.RawTextHelpFormatter
    )

    inputs = parser.add_argument_group("Required Inputs")
    
    inputs.add_argument("-s", "--smiles", required=True,
                        help="Input molecule: SMILES string OR path to .smi file")
    
    inputs.add_argument("-r", "--receptor", required=True,
                        help="Receptor file path (.pdb or .pdbqt)")
    
    inputs.add_argument("--site", action='append', required=True,
                        help="Binding site residues (e.g., --site A:109,110 --site B:210)")
    
    # Optional
    options = parser.add_argument_group("Experiment Options")
    
    options.add_argument("-j", "--job_name", default="untitled",
                         help="Job name for output folder (default: untitled)")
    
    options.add_argument("-c", "--config", default="config/default.yaml",
                         help="Path to config file (default: config/default.yaml)")
    
    options.add_argument("-o", "--output_root", default="results",
                         help="Root output directory (default: results)")

    return parser.parse_args()

# ---------- #
# Run MOZAIC #
# ---------- #

def main():
    args = parse_args()

    # Job name directory
    job_dir = create_job_dir(args.output_root, args.job_name)

    # Log
    log_file = setup_logger(output_dir=job_dir)
    
    start_time = datetime.now() 
    
    logging.info("=" * 50)
    logging.info(f"MOZAIC Started at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logging.info(f"Output Directory  : {job_dir}")
    logging.info(f"Job Name          : {args.job_name}")
    logging.info(f"Log file saved to : {log_file}")
    logging.info("=" * 50 + "\n")
    
    try:
        cfg = load_config(args.config)

        # Config
        shutil.copy(args.config, os.path.join(job_dir, "config.yaml"))
        
        docking_cfg = cfg['docking']
        csa_cfg = cfg['csa']

        # SMILES
        raw_smiles = load_smiles(args.smiles)
        
        mol = Chem.MolFromSmiles(raw_smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES input: {raw_smiles}")
            
        initial_smiles = Chem.MolToSmiles(mol, canonical=True)

        # Receptor
        receptor_path = prepare_receptor(args.receptor, job_dir)

        # Binding center
        binding_center = prepare_center(receptor_path, args.site)
        box_size = docking_cfg['box_size']
        n_cpu    = docking_cfg['n_cpu']
        
        n_rounds = csa_cfg['n_rounds']
        n_bank   = csa_cfg['n_bank']
        n_gen0   = csa_cfg['n_gen0']
        
        logging.info(">>> MOZAIC CONFIGS <<<\n")   
        
        logging.info("-" * 40)
        logging.info(" [1] INPUT DATA")
        logging.info(f"  - Initial Molecule : {initial_smiles}")
        logging.info(f"  - Receptor Path    : {receptor_path}")
        logging.info(f"  - Binding Sites    : {args.site}")
        
        logging.info("-" * 40)
        logging.info(" [2] DOCKING")
        logging.info(f"  - Calculated Center: {binding_center}")
        logging.info(f"  - Box Size         : {box_size}")
        logging.info(f"  - CPU Count        : {n_cpu}")
        
        logging.info("-" * 40)
        logging.info(" [3] CSA")
        logging.info(f"  - Rounds           : {n_rounds}")
        logging.info(f"  - Bank Size        : {n_bank}")
        logging.info(f"  - Gen0 Size        : {n_gen0}")
        logging.info("-" * 40 + "\n")
        
        # Main Logic
        logging.info(">>> MOZAIC START! <<<\n")
        
        final_bank, initial_bank = run_csa(
            job_dir,
            initial_smiles,
            receptor_path,
            binding_center,
            docking_cfg,
            csa_cfg
        )

    except Exception as e:
        logging.error(f"Error: {e}", exc_info=True)
        sys.exit(1)
        
    logging.info(">>> MOZAIC Finished! <<<\n")
    
    end_time = datetime.now()
    duration = end_time - start_time
    
    logging.info("=" * 50 + "\n")
    logging.info(f"Total Runtime: {duration}\n")
    logging.info("=" * 50)

if __name__ == "__main__":
    main()