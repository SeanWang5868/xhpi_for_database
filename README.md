# Characterization of XH-π Interactions in PDB using the Gemmi module

## Introduction

This Python script is used to analyze XH-π interactions based on Plevin system.

## Installation requirements

```bash
pip install gemmi tqdm wget pandas
```

## Usage

### 1. Set `CLIBD_MON` environment variable

**Important:** You need to download the [Monomer library](https://github.com/MonomerLibrary/monomers) and unzip it to your local directory, which is a necessary library for hydrogenating pdb structures.

After unzipping, set the environment variable `CLIBD_MON` in `config.py` to point to the unzipped directory in the code. 

For example, if you unzipped the Monomer library to `/Volumes/Sean/pdb_mirror/monomer/monomers`, you can set the environment variable with the following code:

```bash
os.environ['CLIBD_MON'] = '/Volumes/Sean/pdb_mirror/monomer/monomers'
```

### 2. Run the script

```bash
python xhpi.py
```

### 3. Provide PDB name

```
Enter 4-letter PDB names separated by commas (e.g., 12as,5FJJ,1Gn0):
```

### 4. View results

Once the analysis is complete, the results will be saved to a folder named with the current date. The folder will contain:

## Contact

If you have any questions or need help, please contact Sean Wang (sean.wang@york.ac.uk).
