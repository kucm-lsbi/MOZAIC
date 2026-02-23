import random

from statistics import mean
from .update import distance
    
def make_seed(bank, n_seed):
    
    candidates = [mol for mol in bank if mol['used'] == 'no'] 
    
    seed = random.sample(candidates, n_seed) 
    
    for s in seed: s['used'] = 'yes' 
        
    return seed

def mark_unused(population):

    for item in population:
        item["used"] = "no"
            
    return population

def count_unused(bank):

    return sum(1 for item in bank if item.get("used") == "no")
