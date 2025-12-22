# Fully Automated ECD Workflow for Natural Products

![ORCA](https://img.shields.io/badge/Quantum_Chem-ORCA_6.1.0-blue.svg)
![Python](https://img.shields.io/badge/Python-3.x-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

**An end-to-end automated pipeline for stereochemical determination of flexible natural products (e.g., Lignans) using ORCA and Python.**

## 📖 Overview

This repository contains the workflow described in our paper: *[Insert Paper Title Here]*.

Assigning the absolute configuration of flexible molecules (like lignans) is a challenge due to multiple rotatable bonds. This workflow automates the entire process on **SLURM-managed HPC clusters**, significantly reducing the "wall-clock time" from days to hours (approx. 5-10h for typical lignans).

### Key Features
* **Zero-Intervention:** From 3D conformer search (GOAT/XTB) to TDDFT-ECD calculation in one go.
* **HPC Optimized:** Handles parallel job submission and resource management automatically.
* **Dual-Language Support:** The main script (`orca_double_language.sh`) contains comments in both English and Chinese for easier modification.
* **Automated Visualization:** The Python processor generates publication-ready plots comparing calculated vs. experimental spectra (including Boltzmann weighting and UV-shift correction).

---

## 📂 File Structure

* `orca_double_language.sh`: The master Bash script for SLURM. Controls the entire calculation logic.
* `ecd_processor.py`: Post-processing script for Boltzmann averaging, data normalization, and plotting.
* `split_xyz.py`: (Helper) Splits the ensemble XYZ file into individual conformer files.
* `exp_data.csv`: (User Input) Your experimental ECD data for comparison.

---

## ⚙️ Configuration (Crucial Step)

**Before running the workflow, you MUST modify `orca_double_language.sh` to match your cluster environment.**

Open `orca_double_language.sh` and check the following sections:

### 1. SLURM Header (Lines 15-25)
Modify the `#SBATCH` directives to fit your cluster's policy:
```bash
#SBATCH -p cn-short           # <--- Change to your partition name
#SBATCH -A your_account_name  # <--- Change to your billing account
#SBATCH --qos=your_qos_name   # <--- Change (or remove) based on your QoS

```

### 2. Software Paths (Lines 40-50)

Set the absolute paths to your ORCA and OpenMPI installations:

```bash
# Example:
export ORCA_DIR=/path/to/your/software/orca-6.1.0  # <--- MODIFY THIS
export MPI_HOME=/path/to/your/software/openmpi     # <--- MODIFY THIS

```

### 3. Working Directory (Line 60)

The script attempts to `cd` into a specific directory. Update this to your project folder or remove the line to run in the current submission directory:

```bash
cd /your/project/working/directory  # <--- MODIFY or DELETE THIS

```

---

## 🚀 Usage

### Step 1: Prepare Input Files

1. Place your initial 3D structure (e.g., `structure.xyz`) in the folder.
2. (Optional) Place your experimental data in `exp_data.csv` (Columns: Wavelength, Intensity).

### Step 2: Submit Job

Run the Bash script via SLURM:

```bash
sbatch orca_double_language.sh

```

*The script will automatically perform conformational search, geometry optimization, frequency check, and ECD calculation.*

### Step 3: Post-Processing

Once the job finishes, use the Python script to generate the final comparison plot:

```bash
python3 ecd_processor.py

```

* **Input:** Reads `opt_conf/boltzmann_results.txt` and `exp_data.csv`.
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
* **ORCA 6.1.0** (Required for the specific GOAT/XTB implementation used here).
* **Python 3.x**
* `numpy`
* `matplotlib`



```bash
pip install numpy matplotlib

```

## 📄 Citation

If you use this workflow in your research, please cite:

> [Authors]. "[Title of your Manuscript]". *[Journal Name]*, **Year**, *Vol*, Page. DOI: [Link]

## 📝 License

This project is licensed under the MIT License - see the `LICENSE` file for details.

```
