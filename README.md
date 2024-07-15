
# Characterization of XH-π Interactions in the Protein Data Bank using the Gemmi module

## Introduction

This Python script is used to download CIF files from RCSB PDB, add hydrogen atoms, and analyze XH-π interactions based on Hudson or Plevin methods. The script's functions include:
1. **Download PDB File**: Download the corresponding CIF file based on the 4-letter PDB name entered by the user.
2. **Add Hydrogen Atoms**: Add hydrogen atoms to the PDB file using the `gemmi` tool.
3. **Analyze XH-π Interaction**: Analyze the interaction between hydrogen atoms and aromatic rings using Hudson or Plevin methods.
4. **Handle Missing Monomers**: Record the missing monomer information during the addition of hydrogen atoms.
5. **Result Output**: Save the analysis results as a CSV file.

## Installation requirements

Please make sure the following dependencies are installed in your environment:

- Python 3.x
- `gemmi` library: for processing CIF files and operating structured data.
- `tqdm` library: for displaying progress bars.
- `wget` library: for downloading PDB files.
- `pandas` library: for data processing and result saving.

You can use the following commands to install these dependencies:

```bash
pip install gemmi tqdm wget pandas
```

## Usage

### 1. Set `CLIBD_MON` environment variable

**Important:** You need to download the [Monomer library](https://github.com/MonomerLibrary/monomers) and unzip it to your local directory.

After unzipping, set the environment variable `CLIBD_MON` to point to the unzipped directory in the code. 

For example, if you unzipped the Monomer library to `/path/to/monomers`, you can set the environment variable with the following code:

```bash
os.environ['CLIBD_MON'] = '/Volumes/Sean/pdb_mirror/monomer/monomers'
```
Please replace the path here in the source code.

**Future plans:** We plan to automatically download and manage the Monomer library in future versions to simplify the setup process for users.

### 2. Run the script

Navigate to the directory where the script is located in the terminal and run the following command:

```bash
python script.py
```

### 3. Provide PDB name

After running the script, you will be prompted to enter the 4-letter PDB name. For example:

```
Enter 4-letter PDB names separated by commas (e.g., 12as,5FJJ,1Gn0):
```
Enter your list of PDB names separated by commas, all in uppercase. For example:

```
12AS,5FJJ,1GN0
```

### 4. Select analysis method

Next, you need to select the analysis method:

```
Enter method (Hudson or Plevin):
```

Enter `Hudson` or `Plevin` to select the analysis method. Please note that the entry is case sensitive.

### 5. View results

Once the analysis is complete, the results will be saved to a folder named with the current date. The folder will contain:

- **`missing_monomers.csv`**: records the missing monomer information during the addition of hydrogen atoms.
- **`result_<method>.csv`**: Contains the results of the XH-π interaction analysis, `<method>` is the analysis method you selected (`Hudson` or `Plevin`).

## Example

Here is an example execution flow:

```bash
python XHpi.py
```

Input:

```
Enter PDB names: 1C7H,2J6V,3C6T
Enter method (Hudson or Plevin): Hudson
```

The script will download the CIF files of `1C7H`, `2J6V`, `3C6T`, add hydrogen atoms, and perform XH-π interaction analysis based on the Hudson method. The results will be saved to the `2024-07-15` folder.


## License

The project will be open source after viva.

## Contact Information

If you have any questions or need help, please contact Sean Wang (bql506@york.ac.uk).
