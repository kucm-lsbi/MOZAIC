# MOZAIC 🧩

MOZAIC is a molecule optimization algorithm that combines SMARTS-based reaction-driven fragment growing with the Conformational Space Annealing (CSA) algorithm.

Starting from an initial molecule, MOZAIC performs fragment growing to build an initial bank of candidate molecules. It then applies crossover and mutation operators to increase structural diversity, and updates the bank using the Dcut criterion in CSA. As a result, MOZAIC outputs the globally optimized final bank and corresponding protein–ligand complex PDB structures.

![MOZAIC workflow](./mozaic.png)

---

## 🛠️ Installation & Environment

Clone the repository and create the conda environment:

```bash
git clone https://github.com/kucm-lsbi/MOZAIC.git
cd MOZAIC

conda env create -f environment.yaml
conda activate MOZAIC
```

# 🚀 Usage

## Quick Start

```bash
python main.py \
  -s ./example/IU0.smi \
  -r ./example/7zdn_A.pdb \
  --site A:63,72,73,74,75,76,81,82,101,102,103,107,108,109,110,111 \
  -j 7zdn
```

## Arguments

| Argument | Status | Description |
|---|---|---|
| `-s` | Required | Path to the initial molecule SMILES file (`.smi`). |
| `-r` | Required | Path to the receptor file (`.pdb` or `.pdbqt`). |
| `--site` | Required | Target binding residues formatted as `Chain:ResNum,ResNum,...` (e.g., `A:63,72`). Multiple chains are supported. |
| `-j` | Optional | Job name, used as a prefix for output directories and files. |
| `-c` | Optional | Path to a custom configuration file (default: `./config/default.yaml`). |
| `-o` | Optional | Path to the output directory (default: `./result`). |
| `-h` | Optional | Show the help message and exit. |

# Notes

The output includes the optimized final bank and ranked protein–ligand complex PDB files.

If the receptor is provided in .pdb format, Open Babel must be available in the environment for PDB → PDBQT conversion.
