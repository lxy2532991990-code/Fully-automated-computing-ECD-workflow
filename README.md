# Fully Automated ECD Workflow for Natural Products

![ORCA](https://img.shields.io/badge/ORCA-6.1.0-blue.svg)
![OpenMPI](https://img.shields.io/badge/OpenMPI-4.0.8-orange.svg)
![Python](https://img.shields.io/badge/Python-3.x-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

**An end-to-end automated pipeline for stereochemical determination of flexible natural products (e.g., Lignans) using ORCA and Python.**

## 📖 Overview

This repository contains the workflow described in our paper:

Assigning the absolute configuration of flexible molecules (like lignans) is a challenge due to multiple rotatable bonds. This workflow automates the entire process on **SLURM-managed HPC clusters**, significantly reducing the "wall-clock time" from days to hours (approx. 5-10h for typical lignans).

### Key Features
* **Zero-Intervention:** From 3D conformer search (GOAT/XTB) to TDDFT-ECD calculation in one go.
* **HPC Optimized:** Handles parallel job submission and resource management automatically.
* **Dual-Language Support:** The main script (`orca.sh`) contains comments in both English and Chinese for easier modification.
* **Automated Visualization:** The Python processor generates publication-ready plots comparing calculated vs. experimental spectra (including Boltzmann weighting and UV-shift correction).

---

## 📂 File Structure

* `orca.sh`: The master Bash script for SLURM. Controls the entire calculation logic.
* `ecd_processor.py`: Post-processing script for Boltzmann averaging, data normalization, and plotting.
---

## ⚙️ Configuration (Crucial Step)

**Before running the workflow, you MUST modify `orca.sh` to match your cluster environment.**

Open `orca.sh` and check the following sections:

### 1. SLURM Header (Lines 19-27, 222-230)
Modify the `#SBATCH` directives to fit your cluster's policy:
```bash
#SBATCH -J orca-ecd-calc            # Job name
#SBATCH -p cn-short                 # Partition name (cluster-specific)
#SBATCH -N 1                        # Number of nodes
#SBATCH -n 64                       # Total CPU cores
#SBATCH -o orca-ecd-calc_%j.out     # Standard output file
#SBATCH -e orca-ecd-calc_%j.err     # Standard error file
#SBATCH -A your_account_name        # Replace with your account name
#SBATCH --qos=your_qos_name         # Replace with your Qos (if supported)
#SBATCH --no-requeue                # Do not requeue failed jobs
...
#SBATCH -J orca-conf                # Job name
#SBATCH -p cn-short                 # Partition name (cluster-specific)
#SBATCH -N 1                        # Number of nodes
#SBATCH -n 64                       # Total CPU cores
#SBATCH -o orca-conf_%j.out         # Standard output file
#SBATCH -e orca-conf_%j.err         # Standard error file
#SBATCH -A your_account_name        # Replace with your account name
#SBATCH --qos=your_qos_name         # Replace with your Qos (if supported)
#SBATCH --no-requeue                # Do not requeue failed jobs
```

### 2. Software Paths (Lines 37, 42, 239, 244)

Set the absolute paths to your ORCA and OpenMPI installations:

```bash
export ORCA_DIR=/path/to/your/orca                   # Modify to your ORCA installation path
...
export MPI_HOME=/path/to/your/openmpi                # Modify to your OPENMPI installation path
...
export ORCA_DIR=/path/to/your/orca                   # Modify to your ORCA installation path
...
export MPI_HOME=/path/to/your/openmpi                # Modify to your OPENMPI installation path
```

### 3. Working Directory (Lines 50, 251)

The script attempts to `cd` into a specific directory. Update this to your project folder or remove the line to run in the current submission directory:

```bash
cd /path/to/your/folder                              # Modify to your xyz file folder
...
cd /path/to/your/folder                              # Modify to your xyz file folder
```

### 4. Input Parameter Configuration (Lines 74-75, 149-150, 565)

Configure the following parameters according to your personal situation:

```bash
initial_xyz="your_name_of_structure.xyz"        # Change to your xyz file name
conf_search_base="conf-your_name_of_structure"  # Change to your xyz file name
...
atoms_per_conformer=number_of_atoms             # Number of atoms in molecule
max_conformers_to_process=30                    # Maximum conformers to process, 10 is usually enough for small molecules
...
pool_size=20                                    # Number of concurrently running jobs, adjust according to HPC policy
```

### 5. Calculation Parameter Settings

Set the following parameters based on the characteristics of the compounds and the solvents used:

#### 5.1 Conformational Search (Lines 87-92)

```bash
!GOAT XTB

%pal nprocs 64 end  # Adjust according to your HPC resources
%maxcore 6000       # Adjust according to your HPC resources

* xyzfile 0 1 $initial_xyz
```

#### 5.2 Geometry Optimization (Lines 278-291)

Depending on the properties of your compounds and the available computational resources, you may select different functionals, basis sets, and other parameters.

```bash
! B3LYP OPT FREQ DEF2-TZVP D4 CPCM(Methanol) 
! TightOpt

%pal nprocs 64 end  # Adjust according to your HPC resources
%maxcore 6000       # Adjust according to your HPC resources

%cpcm
end

%freq
   Temp 298.15
end

* xyzfile 0 1 $xyz_file
```

#### 5.3 ECD Calculation (Lines 334-344)

Depending on the properties of your compounds and the available computational resources, you may select different functionals, basis sets, and other parameters.

```bash
!B3LYP DEF2-TZVP CPCM(Methanol)

%pal nprocs 64 end  # Adjust according to your HPC resources
%maxcore 6000       # Adjust according to your HPC resources

%TDDFT
      NROOTS  25
      TDA     FALSE
END

* xyzfile 0 1 opt-conf-${i}.xyz
```
---

## 🚀 Usage

### Step 1: Prepare Input Files

Place your initial 3D structure (e.g., `structure.xyz`) in the folder.

### Step 2: Submit Job

Run the Bash script via SLURM:

```bash
sbatch orca.sh

```

*The script will automatically perform conformational search, geometry optimization, frequency check, and ECD calculation.*

### Step 3: Post-Processing

After the calculation is complete, transfer the `opt_conf` and `ECD_opt` folders to your local computer. Then, place the `ecd_processor.py` and the `exp_data.csv` in the same directory. Finally, modify the `def main()` parameters (Lines 716-733) and run the Python script to obtain the final comparison plot.

```bash
python3 ecd_processor.py

```

* **Input:** Reads `opt_conf/boltzmann_results.txt`, `ECD_opt/ECD-conf-*.out` and `exp_data.csv`.
* **Output:** Generates `comparison_plot.png` and `calculated_spectrum.csv`.

---

## 📊 Output Examples

The workflow produces organized outputs:

1. **`opt_conf/`**: Contains optimized geometries (`.xyz`) and Boltzmann distribution data (`boltzmann_results.txt`).
2. **`ECD_opt/`**: Contains raw ORCA output files for every conformer.
3. **Visualizations**:
* **Boltzmann Averaged Spectrum**: Adjusted for conformer population.
* **Comparison Plot**: Overlays Experimental (Black) vs. Calculated (Red/Blue for enantiomers).



---

## 🛠️ Requirements

* **Linux HPC Cluster** with SLURM workload manager.
* **ORCA 6.1.0**
* **OpenMPI 4.0.8**
* **Python 3.x**
* `numpy`
* `matplotlib`



```bash
python -m pip install numpy matplotlib

```

## 📄 Citation

If you use this workflow in your research, please cite:

## 📝 License

This project is licensed under the MIT License - see the `LICENSE` file for details.
