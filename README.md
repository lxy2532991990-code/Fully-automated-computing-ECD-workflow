### Fully-automated-computing-ECD-workflow

Automated workflow for computing ECD spectra with ORCA. Features an end-to-end pipeline from conformational search to Boltzmann-weighted spectral generation and comparison with experimental data. Includes tools for parallel HPC job management and automated data processing. Designed to save time and reduce errors for computational chemists.

## 🚀 Features

* **End-to-End Automation**: Handles the complete lifecycle from 3D conformer search (GOAT/XTB) to final ECD spectrum generation.
* **HPC Integration**: Optimized `SLURM` submission script (`orca.sh`) that manages parallel execution of geometry optimizations and TDDFT calculations for multiple conformers.
* **Smart Post-Processing**: Python-based processor (`ecd_processor.py`) that:
* Calculates Boltzmann weights based on Gibbs free energy.
* Applies Gaussian broadening to discrete transitions.
* Performs UV-shift correction and intensity scaling.
* Supports automatic spectrum inversion for enantiomer determination.


* **Publication-Ready Visualization**: Automatically generates high-quality comparison plots (Calculated vs. Experimental) using `matplotlib`.

## 🛠️ Requirements

### Software

* **ORCA 6.1.0** (or later)
* **OpenMPI** (compatible with your ORCA version)
* **Python 3.x**

### Python Dependencies

Install the required libraries via pip:

```bash
pip install numpy matplotlib

```

## 📂 File Structure

* `orca.sh`: Main Bash script for SLURM workload manager. Handles environment setup, conformational search, geometry optimization, and ECD calculation steps.
* `ecd_processor.py`: Post-processing tool for Boltzmann weighting, spectral fitting, and plotting.
* `exp_data.csv`: (User provided) Experimental ECD data file containing two columns: Wavelength (nm) and Intensity (mdeg or Δε).

## 📖 Usage

### Step 1: HPC Calculation (`orca.sh`)

1. Modify the header of `orca.sh` to match your cluster's partition, account, and resource limits.
2. **Important:** Update the `ORCA_DIR` and `MPI_HOME` paths in the script to point to your local installation:
```bash
export ORCA_DIR=/path/to/your/orca
export MPI_HOME=/path/to/your/openmpi

```


3. Submit the job:
```bash
sbatch orca.sh

```


*The script will automatically perform conformational search, optimize geometries for all conformers, and calculate ECD/UV spectra.*

### Step 2: Post-Processing (`ecd_processor.py`)

Once the calculations are complete, use the Python script to generate the final spectra.

1. Prepare your experimental data in `exp_data.csv`.
2. Run the processor:
```bash
python3 ecd_processor.py

```


3. **Parameter Tuning:** You can adjust the spectral parameters inside the script or (future version) via command line arguments:
* `sigma`: Standard deviation for Gaussian broadening (default: 0.3 eV).
* `shift`: UV-shift correction in nm (default: 0).
* `scale_factor`: Scaling factor for intensity matching.



## 📊 Outputs

The workflow generates the following key files:

* **`boltzmann_results.txt`**: Detailed list of conformer energies and their calculated Boltzmann populations.
* **`calculated_spectrum.csv`**: Raw and scaled calculated ECD data points.
* **`comparison_plot.png`**: A visual overlay of the experimental spectrum against the Boltzmann-averaged calculated spectrum (including inverted configuration for comparison).

## 📝 Citation

If you use this workflow in your research, please cite:

> [Your Name], et al. "Title of your Manuscript". *Journal Name*, Year, Volume, Pages. (DOI)

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.
