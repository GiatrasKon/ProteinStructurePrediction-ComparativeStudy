
# Protein Structure Prediction: Rosetta Ab Initio vs. Linear Algebra Algorithm

This repository contains the code, data, and documentation for a project comparing the performance of two computational methods for protein structure prediction: **[Rosetta Ab Initio](https://docs.rosettacommons.org/docs/latest/application_documentation/structure_prediction/abinitio)** and a **Linear Algebra Algorithm**, described in the ["Molecular conformation search by distance matrix perturbations" by Emiris I. et. al](https://link.springer.com/article/10.1007/s10910-004-1466-4). These methods are critical for understanding molecular mechanisms, advancing drug design, and exploring diseases at a molecular level.

This analysis was performed as part of the final project for the "Algorithms in Structural Bioinformatics" graduate course of the MSc Data Science & Information Technologies Master's programme (Bioinformatics - Biomedical Data Science Specialization) of the Department of Informatics and Telecommunications department of the National and Kapodistrian University of Athens (NKUA), under the supervision of professors Ioannis Emiris and Evangelia Chrysina, in the academic year 2022-2023.

---

## Contributors
- [Konstantinos Giatras](https://github.com/GiatrasKon)
- [Spyros Alvanakis-Apostolou](https://github.com/SpyrosAlvanakis)

---

## Workflow and Tools

1. **Rosetta Ab Initio**
    - Predicts tertiary protein structures from amino acid sequences.
    - Relies on Monte Carlo simulations and energy minimization.
    - Input: FASTA sequence, secondary structure predictions, and fragment libraries.
    - Output: Predicted protein structures as PDB files.
2. **Linear Algebra Algorithm**
    - Processes Cayley-Menger matrices derived from NMR distance data.
    - Perturbs matrices iteratively using Singular Value Decomposition (SVD).
    - Output: Symmetric matrices for extraction of 3D coordinates.
3. Supplementary Tools
    - Python Scripts: For RMSD calculations and handling predicted structures.
    - MATLAB: For executing the Linear Algebra algorithm.
---

## Results Summary

**Rosetta Ab Initio**:
    - More effective for small-to-medium-sized proteins.
    - Accuracy improves with the number of predicted conformations (10, 100, 1000 tested).
**Linear Algebra Algorithm**:
    - Limited to small proteins due to high computational complexity (O(n³) for SVD).
    - Did not yield accurate results for larger proteins due to constraints in matrix rank reduction.

---

## Installation & Usage

### Cloning the Repository

```sh
git clone https://github.com/GiatrasKon/ProteinStructurePrediction-ComparativeStudy.git
```

### Package Dependencies

Ensure you have the following packages installed:

- pandas
- numpy
- scikit-learn
- biopython
- sklearn
- matplotlib
- seaborn

Install dependencies using:

```sh
pip install pandas numpy biopython scikit-learn matplotlib seaborn 
```

---

### Repository Structure

- `documents/`: Project documentation and references
    - `AiSB - Final Project Report - Alvanakis-Giatras.pdf`: Final project report
    - `AiSB - Final Presentation - Alvanakis-Giatras.pdf`: Final project presentation
    - `EmirNiki.pdf`: Reference paper for the Linear Algebra Algorithm
- `Linear Algebra Algorithm/`: Implementation of the numerical rank-reduction approach
    - `data/`: Input data for Matlab processing (bzip2 compressed, need to be extracted running `tar -xvjf data.tar.bz2`)
        - `*.pdb`: Protein structure files (e.g., 1adx.pdb, 1anp.pdb)
        - `*.m`: Cayley-Menger matrices and bounding matrices
        - `*.csv`: Ouput Cayley-Menger matrices and bounding matrices
    - `scripts/`: Matlab scripts implementing the Linear Algebra Algorithm
        - `inout/`: Example input/output files
        - `matlab/`: Core algorithm implementation
            - `mconf.m`: Main script for Cayley-Menger matrix rank reduction
            - `svred.m`: Singular value decomposition optimization
            - `embed.m`: Embedding distance matrices into 3D space
            - Additional utility functions
    - `notebooks/`
        - `LAA_Supplementary_Code.ipynb`: Supplementary Python code for processing Matlab outputs
- `Rosetta Ab Initio/`: Implementation of Rosetta-based protein structure prediction
    - `data/`: Input and output data
        - `<PDB_ID>_pdbs/`: Directory for each protein (e.g., 2jyv_pdbs, 2mod_pdbs)
            - `10/`: Contains PDBs for 10 predicted structures
            - `100/`: Contains PDBs for 100 predicted structures
            - `1000/`: Contains PDBs for 1000 predicted structures
            - `ensemble_states/`: Isolated states of each protein ensemble
            - `aat000_03_05.200_v1_3`: 3-mer fragment library file for structure prediction
            - `aat000_09_05.200_v1_3`: 9-mer fragment library file for structure prediction
            - `*.fasta`: Protein sequence files
            - `*.out`: Silent files with predicted structures
            - `flags`: Configuration files for Rosetta execution
            - `*.ss2`: Secondary structure prediction files
        - `top_*.csv`: Results with top predictions (e.g., RMSD values)
    - `notebooks/`:            
        - `rename_silent.ipynb`: Supplementary scripts for Rosetta processing
        - `RAI_Supplementary_Code.ipynb`: Supplementary Python code for RMSD calculations

---

### Usage

This section provides step-by-step instructions to recreate the analysis for both the **Linear Algebra Algorithm** and **Rosetta Ab Initio** methods. Follow these steps in the specified order to ensure accurate reproduction of the results.


#### 1. General Setup

Before beginning, ensure you have the following installed:
- **MATLAB** (for running the Linear Algebra Algorithm, with the required toolboxes e.g. Linear Algebra Toolbox for SVD)
- **Python 3.x** with the necessary packages (refer to the Python package dependencies section)
- **PyMOL** (for visualizing PDB structures)
- **Rosetta software suite** (specifically the AbinitioRelax tool, properly licensed and installed)

##### Software-Specific Notes:
- **MATLAB**: Ensure it can run `.m` files in the specified working directory.
- **Rosetta**: Obtain the software from the official [RosettaCommons](https://www.rosettacommons.org/software/) website and acquire the necessary license. The `flags` files and protein-specific input files (e.g., `.fasta`, `.ss2`, 3-mer and 9-mer fragment files) are essential.

##### Data Requirements:
- Protein structures should be downloaded from the Protein Data Bank (PDB) in `.pdb` format. For this project, the following PDB entries were used:
  - 2JYV, 2M0D, 2MPC, 2MS8, 6MWM for Rosetta Ab Initio.
  - 1ADX and 1ANP for the Linear Algebra Algorithm.


#### 2. Rosetta Ab Initio

##### Input Preparation:
1. **Download Protein Structures**:
   - Obtain PDB files for your proteins from the [PDB website](https://www.rcsb.org/).
   - Save each protein's `.pdb` file in a separate directory (e.g., `Rosetta Ab Initio/data/<PDB_ID>_pdbs/`).

2. **Generate Fragment Files**:
   - Use the [Robetta Fragment Server](http://robetta.bakerlab.org/fragmentsubmit.jsp) to generate the 3-mer (`aat000_03_05.200_v1_3`) and 9-mer (`aat000_09_05.200_v1_3`) fragment files, as well as the secondary structure prediction file (`t000_.psipred_ss2`).
   - Save these files in the corresponding protein directory.

3. **Prepare a Flags File**:
   - Create a `flags` file for each protein, specifying paths to the `.pdb`, `.fasta`, `.ss2`, and fragment files. Below is an example template:
```sh
-database /path/to/rosetta/main/database
-in:file:native ./<PDB_ID>.pdb
-in:file:fasta ./<PDB_ID>.fasta
-in:file:frag3 ./aat000_03_05.200_v1_3
-in:file:frag9 ./aat000_09_05.200_v1_3
-psipred_ss2 ./t000_.psipred_ss2
-nstruct 1000
-abinitio:relax
-use_filters true
-abinitio::increase_cycles 10
-abinitio::rg_reweight 0.5
-abinitio::rsd_wt_helix 0.5
-abinitio::rsd_wt_loop 0.5
-relax::fast
-out:file:silent ./<PDB_ID>_silent_1000_tries.out
```

##### Running Rosetta:
1. Navigate to the Rosetta directory:
```sh
cd /path/to/rosetta/main/source/bin
```
2. Execute the `AbinitioRelax` tool using the flags file:
```sh
./AbinitioRelax.default.linuxgccrelease @/path/to/<PDB_ID>_pdbs/flags
```

##### Post-Processing:
1. Use the `rename_silent.ipynb` script to extract PDB files from the silent files and rename them for better organization.
2. Run the `RAI_Supplementary_Code.ipynb` notebook to calculate RMSD values for the predicted structures and identify the best conformations.


#### 3. Linear Algebra Algorithm

##### Input Preparation:
1. **Download Protein Structures**:
   - Obtain PDB files for proteins (e.g., 1ADX, 1ANP).
   - Save these in the `Linear Algebra Algorithm/data/` directory.

2. **Generate Cayley-Menger Matrices**:
   - Use the `LAA_Supplementary_Code.ipynb` notebook to:
     - Extract backbone coordinates from the PDB files.
     - Add noise (up to 2%) to generate perturbation matrices simulating NMR data.
     - Construct Cayley-Menger matrices and save them as `.m` files.

##### Running the Algorithm:
1. Open MATLAB and set the working directory to `Linear Algebra Algorithm/scripts/matlab/`.
2. Ensure all required `.m` files and input matrices are present in the directory.
3. Run the main script:
   - `mconf('cm_mtrx_<PDB_ID>.m', 1e-8, 1);`
4. Save the output Cayley-Menger matrix and verify its rank using:
   - `rank(ans);`

##### Post-Processing:
1. Convert the Cayley-Menger matrix back into a distance matrix and 3D coordinates using `LAA_Supplementary_Code.ipynb`.
2. Calculate the RMSD between the predicted coordinates and ground truth using the same notebook.


#### 4. Results Analysis
- For both methods, predicted structures and RMSD calculations will be stored in the respective `data/<PDB_ID>_pdbs/` or `Linear Algebra Algorithm/data/` directories.
- Use PyMOL to visualize the best conformations and compare them with the ground truth.

---