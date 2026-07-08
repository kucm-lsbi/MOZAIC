import copy

from ..scoring.scorer import get_scores

from .crossover import run_crossover
from .mutation import run_mutation

def make_child(initial_bank, current_bank, seed, initial_smiles, receptor_path, binding_center, docking_cfg, csa_cfg):

    child_cfg = csa_cfg.get('child', {})

    initial_bank_copy = copy.deepcopy(initial_bank)
    current_bank_copy = copy.deepcopy(current_bank)
    seed_copy = copy.deepcopy(seed)

    # Crossover: len(seed) * crossover_per_seed
    crossover_population, backup_crossover = run_crossover(
        seed_copy,
        initial_bank_copy,
        current_bank_copy,
        initial_smiles,
        child_cfg
    )

    # Mutation: fills up to len(seed) * population_per_seed
    child, backup_child = run_mutation(
        crossover_population,
        backup_crossover,
        seed_copy,
        initial_smiles,
        child_cfg
    )

    # Scoring
    scored_child = get_scores(
        child, 
        backup_child, 
        receptor_path, 
        binding_center, 
        docking_cfg
    )
    
    return scored_child
