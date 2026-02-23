import json
import logging
from pathlib import Path
from typing import Dict, List, Any

BASE_DIR = Path(__file__).resolve().parents[2] # ProjectRoot
DATA_DIR = BASE_DIR / "data"

FRAG_DIR = DATA_DIR / "fragments_position"
FG_JSON  = DATA_DIR / "reactions" / "functional_groups.json"
RXN_JSON = DATA_DIR / "reactions" / "reactions.json"

class ReactionLibrary:
    
    _instance = None
    
    fragments: Dict[str, List[str]]
    
    functional_groups: Dict[str, str]
    rxn_rules: Dict[str, Any]

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(ReactionLibrary, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if self._initialized:
            return
        self._load_data()
        self._initialized = True

    def _load_data(self):
        logging.info("Loading Reaction Library and Fragments...")
        
        # Load Fragments
        self.fragments = {}
        if FRAG_DIR.exists():
            for smi_path in FRAG_DIR.glob("*.smi"):
                try:
                    with smi_path.open("r", encoding='utf-8') as f:
                        smiles_list = [
                            line.split()[0].strip() 
                            for line in f if line.strip()
                        ]
                        if smiles_list:
                            self.fragments[smi_path.stem] = smiles_list
                except Exception as e:
                    logging.error(f"Failed to load fragment {smi_path.name}: {e}")
                    
            if not self.fragments:
                logging.error(f"No fragment files found in {FRAG_DIR}")
                
        else:
            logging.error(f"Fragment directory not found: {FRAG_DIR}")

        # Load Functional Groups
        self.functional_groups = {}
        if FG_JSON.exists():
            try:
                with FG_JSON.open("r", encoding='utf-8') as f:
                    self.functional_groups = json.load(f)
                    
            except Exception as e:
                logging.error(f"Failed to load Functional Groups: {e}")
        else:
            logging.error(f"Functional Group file not found: {FG_JSON}")

        # Load Reaction Rules
        self.rxn_rules = {}
        if RXN_JSON.exists():
            try:
                with RXN_JSON.open("r", encoding='utf-8') as f:
                    self.rxn_rules = json.load(f)
            except Exception as e:
                logging.error(f"Failed to load Reaction Rules: {e}")
        else:
            logging.error(f"Reaction Library file not found: {RXN_JSON}")
            
        logging.info(f"Library Loaded: {len(self.fragments)} fragment sets, {len(self.rxn_rules)} reactions.")


LIBRARY = ReactionLibrary()

SMILES_FRAGMENTS  = LIBRARY.fragments
FUNCTIONAL_GROUPS = LIBRARY.functional_groups
RXN_RULES         = LIBRARY.rxn_rules