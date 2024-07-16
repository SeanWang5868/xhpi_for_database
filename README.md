# Characterization of XH-π Interactions in the Protein Data Bank using the Gemmi module

## Introduction

This Python script is used to analyze XH-π interactions based on Hudson or Plevin system.
1. **Download PDB File**: Download the corresponding CIF file based on the 4-letter PDB name.
2. **Add Hydrogen Atoms**: Add hydrogen atoms to the PDB file using `gemmi`, record the missing monomers during the addition.
3. **Analyze XH-π Interaction**: Analyze the XH-π Interaction using Hudson or Plevin system.
4. **Result Output**: Save the analysis results as a CSV file.

## Installation requirements

```bash
pip install gemmi tqdm wget pandas
```

## Usage

### 1. Set `CLIBD_MON` environment variable

**Important:** You need to download the [Monomer library](https://github.com/MonomerLibrary/monomers) and unzip it to your local directory.

After unzipping, set the environment variable `CLIBD_MON` in `config.py` to point to the unzipped directory in the code. 

For example, if you unzipped the Monomer library to `/path/to/monomers`, you can set the environment variable with the following code:

```bash
os.environ['CLIBD_MON'] = '/Volumes/Sean/pdb_mirror/monomer/monomers'
```

### 2. Run the script

```bash
python script.py
```

### 3. Provide PDB name and  select the detection system

```
Enter 4-letter PDB names separated by commas (e.g., 12as,5FJJ,1Gn0):
```

```
Select the detection system (Hudson or Plevin):
```

### 4. View results

Once the analysis is complete, the results will be saved to a folder named with the current date. The folder will contain:

- **`ABCD.cif.gz`**: the cif files that just download
- **`missing_monomers.csv`**: records the missing monomer information during the addition of hydrogen atoms.
- **`result_<system>.csv`**: Contains the results of the XH-π interaction analysis, `<system>` is the analysis system selected (`Hudson` or `Plevin`).

## Example

```
Output directory created: ./2024-07-17
Enter 4-letter PDB names separated by commas (e.g., ABCD,EFGH,IJKL):
5fjj,12as
PDB to be processed: 5FJJ, 12AS
100% [............................................................................] 761830 / 761830 0/2 [00:00<?, ?it/s]
Downloaded and saved 5FJJ to ./2024-07-17/5FJJ.cif.gz
100% [............................................................................] 143540 / 143540:00<00:00,  1.35it/s]
Downloaded and saved 12AS to ./2024-07-17/12AS.cif.gz
Downloading PDB files: 100%|██████████████████████████████████████████████████████████████| 2/2 [00:01<00:00,  1.57it/s]
Select the detection system (Hudson or Plevin): Hudson
Processing files:   0%|                                                                           | 0/2 [00:00<?, ?it/s]
Successfully added hydrogen atoms for file ./2024-07-17/5FJJ.cif.gz.
Successfully saved the modified structure to ./2024-07-17/5FJJ.cif.gz.
Processing files:  50%|█████████████████████████████████▌                                 | 1/2 [00:00<00:00,  1.24it/s]
Successfully added hydrogen atoms for file ./2024-07-17/12AS.cif.gz.
Successfully saved the modified structure to ./2024-07-17/12AS.cif.gz.
Processing files: 100%|███████████████████████████████████████████████████████████████████| 2/2 [00:09<00:00,  4.70s/it]
Detect 642 XH-π interactions totally
The result has been saved to ./2024-07-17/result_hudson.csv

```

## License

The project will be open source after viva.

## Contact Information

If you have any questions or need help, please contact Sean Wang (bql506@york.ac.uk).
