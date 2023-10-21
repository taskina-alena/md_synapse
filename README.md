
# MD Simulation of Synaptic Vesicles

This project involves Molecular Dynamics (MD) simulation of synaptic vesicles. The steps to run and analyze the simulations are outlined below.

## 1. Preparation: Create Input Files

Use the `make_simulation.py` script to generate input files required for LAMMPS simulations.

```bash
python make_simulation.py
```

This will create necessary input files in the following directories:
- `Input/Scripts`
- `Input/Configuration`

## 2. Running MD Simulations with LAMMPS

Navigate to the LAMMPS executable directory and run the simulation using:

```bash
../path/to/lmp_serial -in Input/Scripts/*.in
```

**Note:** 
- The simulation results, including movies and logs, will be saved in the `Results` folder.
- You can view the simulation movies using [Ovito](https://www.ovito.org/).

## 3. Analysis

To analyze the movies, find clusters, compute their sizes, and determine their circularity, use the `run_analysis.py` script:

```bash
python Analysis/run_analysis.py
```

The results will be saved in the `Results` folder. Files containing the analysis start with "Oligo...".

## 4. Visualization

You can visualize the results and perform further analysis using the provided Jupyter notebook:

```bash
jupyter notebook Analysis/analysis_plots.ipynb
```
