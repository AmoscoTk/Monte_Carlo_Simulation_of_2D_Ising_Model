# Monte_Carlo_Simulation_of_2D_Ising_Model

### Simulating phase transitions in 2D spin systems via Monte Carlo (Metropolis) — coded in MATLAB

#### 📖 Overview

This project implements a Monte Carlo simulation of the two-dimensional (2D) Ising model using the Metropolis algorithm. The goal is to investigate phase transitions (ferromagnetism → paramagnetism) and thermodynamic quantities such as magnetization, energy, specific heat, susceptibility, etc.

**Key motivations:**
	•	Study how spins on a lattice interact and align at different temperatures
	•	Observe critical behavior (e.g., critical temperature, fluctuations)
	•	Visualize spin configurations and their evolution
	•	Provide a flexible, extensible MATLAB toolkit for educational or research purposes

#### 🧠 Features
	•	Implementation of the Metropolis algorithm for 2D Ising lattices
	•	Ability to tune lattice size, temperature, number of Monte Carlo steps
	•	Thermalization routines to discard initial transient states
	•	Data collection for thermodynamic quantities:
  - Magnetization & magnetization fluctuations
  - Energy & energy fluctuations
  - Heat capacity, magnetic susceptibility
	•	Visualization of spin configurations and time evolution
	•	Modular script structure (in /Scripts, /Visualization, etc.)
	•	Report and analysis files to demonstrate results from simulations

#### 📂 Project Structure

/
├── Scripts/                   # Main simulation scripts  
├── Visualization/             # Plotting & visualization modules  
├── Thermalization time/       # Supporting scripts for equilibration  
├── Report/                    # Analyses, plots, report material  
├── LICENSE (BSD-3-Clause)  
└── README.md                  # (this file)

	•	Scripts/ — simulation entry point(s) (e.g. main runner, parameter sweeps)
	•	Visualization/ — routines for plotting spin lattice, energy / magnetization over time, etc.
	•	Thermalization time/ — tests / routines to estimate how many Monte Carlo steps to discard before measurements
	•	Report/ — sample results, graphs, discussion

#### 🚀 Usage

Here’s how to use the simulation:
	1.	Clone the repo:

git clone https://github.com/AmoscoTk/Monte_Carlo_Simulation_of_2D_Ising_Model.git
cd Monte_Carlo_Simulation_of_2D_Ising_Model


	2.	Open MATLAB and navigate to the project directory.
	3.	Configure simulation parameters in one of the scripts (e.g. lattice size L, temperature range T, number of Monte Carlo steps nSteps, etc.).
	4.	Run the simulation script(s) (e.g. run_simulation.m or similar).
	5.	After simulation, use the visualization scripts to plot:
	•	Spin configuration snapshots
	•	Magnetization / energy vs time or vs temperature
	•	Derived quantities like heat capacity or susceptibility
	6.	Refer to the Report/ directory for example plots and analysis.

#### 🧩 Configuration Parameters (Examples)

Parameter	Meaning / Use	Typical Value / Notes
L	Linear dimension of 2D lattice (L×L spins)	e.g. 20, 50, 100
T	Temperature (in units where k_B = 1)	Sweep over e.g. 1.0 to 4.0
nSteps	Number of Monte Carlo iterations	e.g. 1e5, 5e5
nThermalSteps	Number of steps to discard (thermalization)	e.g. 1e4
measureInterval	Interval of steps between measurements	e.g. every 10 steps

You can adjust these in your simulation scripts to explore different regimes (small vs large systems, near critical T, etc.).

#### 🧪 Example Results & Insights
	•	At low T, spins tend to align (spontaneous magnetization)
	•	As T increases past a critical point T_c ≈ 2.269 (for infinite 2D Ising), magnetization drops to zero
	•	Near T_c, fluctuations in energy and magnetization become large — observable peaks in heat capacity and susceptibility
	•	Finite size effects are present: in finite lattices, transitions are smoothed and critical peaks broaden

You’ll find example plots and discussions in the Report/ folder to help understand these behaviors.

#### 🎯 Potential Extensions

You might consider extending this project with:
	•	Wolff or Swendsen–Wang cluster algorithms (often more efficient near criticality)
	•	Simulations of the 3D Ising model
	•	External magnetic field h ≠ 0
	•	Parallel implementation (GPU or multi-threaded)
	•	Finite-size scaling analysis to estimate critical exponents
	•	Monte Carlo over disorder (e.g. random field Ising model)

#### 📚 References & Theory

Some useful background reading:
	1.	K. Huang, Statistical Mechanics
	2.	Newman & Barkema, Monte Carlo Methods in Statistical Physics
	3.	Classic articles on the 2D Ising solution (Onsager)
	4.	Tutorials on the Metropolis algorithm

#### 📝 License & Contribution

This project is licensed under the MIT License. Any contributions, bug-fixes, or improvements are welcome.