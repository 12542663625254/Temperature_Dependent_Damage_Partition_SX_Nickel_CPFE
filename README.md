# Temperature-Dependent Damage Partition in Single-Crystal Nickel (CPFEM)

A dislocation density-based **crystal plasticity finite element (CPFEM)** model for the creep–fatigue response of single-crystal nickel superalloy (DD6). The model predicts creep–fatigue life with validated **temperature dependence** and produces a **damage-partition map** that separates creep- and fatigue-dominated damage across a wide range of loading configurations. The framework is implemented in **MOOSE**(https://mooseframework.inl.gov/). 


---

## Table of contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Dependencies](#dependencies)
- [Run sequence](#run-sequence)
  - [1. Build the MOOSE application](#1-build-the-moose-application)
  - [2. Material, crystallography and orientation files](#2-material-crystallography-and-orientation-files)
  - [3. Select the loading case and boundary conditions](#3-select-the-loading-case-and-boundary-conditions)
  - [4. Run the simulation](#4-run-the-simulation)
  - [5. Postprocessing](#5-postprocessing)
- [Input files](#input-files)
- [Selecting the loading case](#selecting-the-loading-case)
- [Postprocessing data — naming convention](#postprocessing-data--naming-convention)
- [Postprocessing sub-folders](#postprocessing-sub-folders)
- [Damage calculation](#damage-calculation)

---

## Overview

The custom material implements a dislocation density-based crystal plasticity constitutive model: a power-law flow rule on the 12 FCC slip systems, dislocation-density evolution, an Armstrong–Frederick back-stress with dynamic and Arrhenius static recovery, a fully implicit local Newton–Raphson stress update, and a consistent algorithmic tangent. The simulation domain is a 1 µm cubic representative volume element (RVE) meshed with HEX8 elements.

---

## Repository structure

| Folder / file | Description |
|---|---|
| `src/materials` | MOOSE C++ source of the custom crystal plasticity material model. |
| `include/materials` | Header files (`.h`) for the custom material model. |
| `Input file for MOOSE` | MOOSE input files (`.i`) and the properties, slip system, crystallography input files (`.in`). |
| `Different Loading condition` | Input-file variants for the loading cases (pure fatigue, tensile/compressive hold, balanced 30/30, R-ratio and temperature variants). |
| `Matlab Code and Data for postprocessing of the results` | MATLAB postprocessing scripts and the CSV data used to generate the figures. |
| `README.md` | This file. |

**Source and header files**

| File | Description |
|---|---|
| `RVE_Paper_SXNickel.C` / `.h` | The custom MOOSE material (UMAT-equivalent stress update). Implements the flow rule on the 12 FCC slip systems, forest-dislocation-density evolution, Armstrong–Frederick back-stress with dynamic and Arrhenius static recovery, the implicit local Newton–Raphson update, and the algorithmic tangent. Registered in the app as object type `RVE_Paper_SXNickel`. |

---

## Dependencies

| Component | Purpose |
|---|---|
| [MOOSE](https://mooseframework.inl.gov/) | Finite element framework; the model is a custom MOOSE material. |
| [PETSc](https://petsc.org/) | Linear algebra / solvers (bundled with MOOSE). |
| rho-CP (Patra et al., 2023) | Dislocation density-based CP framework the material builds on. |
| MATLAB | Postprocessing and figure generation. |

---

## Run sequence

The run sequence consists of the following steps:

1. Build the MOOSE application with the custom material.
2. Set up the material, crystallography and orientation input files.
3. Select the loading case and boundary conditions in the input file.
4. Run the simulation.
5. Postprocess the results in MATLAB.

### 1. Build the MOOSE application

Source: `src/materials/RVE_Paper_SXNickel.C` and header `include/materials/RVE_Paper_SXNickel.h`.

Copy the source and header into your MOOSE application's source tree (`src/materials` and `include/materials`), then compile:

```bash
make -j <number_of_processors>
```

Prerequisites: a working MOOSE installation with PETSc. The material registers itself in the application under the object type `RVE_Paper_SXNickel`.

### 2. Material, crystallography and orientation files

These files are read by the `[Materials]` block of the input file. Set them up before running.

| File | Description |
|---|---|
| `fcc_props.in` | The 30 material parameters, one per line, in the order read by `RVE_Paper_SXNickel.C` (elastic constants C11, C12, C44; shear modulus; Burgers vector; hardening/recovery parameters; dislocation densities; strain-rate sensitivity `m`; etc.). **Edit this file to change any material parameter.** |
| `fcc_slip_sys.in` | The 12 FCC slip systems; each row = slip-plane normal (3 integers) followed by slip direction (3 integers). Not normally changed for FCC. |
| `orientations.in` | The single-crystal orientation. Line 1 = number of orientations; following line = Euler angles (`phi1 Phi phi2`). Change to study a different orientation. |

### 3. Select the loading case and boundary conditions

Code: `Input file for MOOSE/fatigue_test.i`

Required changes:

- **Mesh** — 1 µm cube with a 10×10×10 HEX8 mesh (`[Mesh]` block). Change `nx, ny, nz` for a mesh-convergence study.
- **Loading waveform** — in the `[Functions]` block, uncomment the `top_push` block for the desired case and comment out the others.
- **Strain amplitude** — set by the y-values of the waveform (±0.012 = 1.2% strain).
- **R-ratio / mean strain** — set by shifting the waveform y-values (e.g. tension-biased for R > −1).
- **Temperature** — set by `temp` in the `[Materials]` block (1033 K shown; use the value for 760 °C or 980 °C as required).
- **Boundary conditions** — displacement-controlled: the loaded face is driven by `top_push` in z; opposite faces are constrained (`disp = 0`) to prevent rigid-body motion while allowing Poisson contraction.
- **Time stepping** — the `dts` function and `[Executioner]` block control the time step; small steps are used at the onset of hold periods for the stress-relaxation kinetics.

### 4. Run the simulation

```bash
mpirun -np <number_of_processors>  <your-app>-opt  -i "Input file for MOOSE/fatigue_test.i"
```

Replace `<your-app>-opt` with the name of your compiled MOOSE application. The simulation writes the stabilized-cycle stress, strain, dislocation density, back-stress and other fields to CSV/Exodus output as defined in the `[Outputs]` block.

### 5. Postprocessing

Folder: `Matlab Code and Data for postprocessing of the results`

- **Collect the CSV output** — name the file by the [naming convention](#postprocessing-data--naming-convention).
- **Run the MATLAB script** — it reads the CPFEM CSV output, extracts `strain_zz` (column 33) and `stress_zz` (column 36), and plots the stress-strain hysteresis response.
- **Output** — the script produces the stress-strain figure (with a zoomed inset of the hardening / hold-relaxation region). These correspond to the stress-strain and hysteresis figures in the paper.

---

## Input files

| File | Description |
|---|---|
| `fatigue_test.i` | The main MOOSE input file. Defines the RVE (10×10×10 HEX8), boundary conditions, loading waveform, material block, solver settings and output. The loading case is selected inside this file. |
| `fcc_props.in` | The 30 material parameters (see [step 2](#2-material-crystallography-and-orientation-files)). |
| `fcc_slip_sys.in` | The 12 FCC octahedral slip systems (normal + direction per row). |
| `orientations.in` | Crystal orientation (Euler angles). `0 0 0` = cube-oriented single crystal. |

---

## Selecting the loading case

The `[Functions]` block of `fatigue_test.i` contains the strain-vs-time waveform (the `top_push` function) for every loading case, one after another. **All but one are commented out.** To run a given case, uncomment its `top_push` block and comment out the others.

| Case | Description |
|---|---|
| Pure fatigue (0/0) | Continuous triangular waveform; provided at 1.0%, 1.2% and 1.5% strain amplitude. |
| Tensile hold (60/0) | 60 s dwell at the tensile peak each cycle; also 30/120/240/300 s tensile holds. |
| Compressive hold (0/60) | 60 s dwell at the compressive peak each cycle. |
| Balanced hold (30/30) | 30 s dwell at both peaks. |
| Stress-amplitude cyclic tests | Fully reversed cycling for the stress-amplitude-vs-cycle plots (1.0% and 1.5%). |

---

## Postprocessing data — naming convention

CSV data files are named by the pattern:

```text
<hold_config>_CPFEM_upto_<strain>_percent_<cycle>_cycle.csv
```

| Field | Meaning |
|---|---|
| `hold_config` | Dwell configuration: `30_30` (balanced), `60_0` (tensile hold), `0_60` (compressive hold), or `0_0` / absent for pure fatigue. |
| `CPFEM` | Model identifier (dislocation density-based crystal plasticity, implicit). |
| `strain` | Strain amplitude, underscore for the decimal point, e.g. `1_2_percent` = 1.2%. |
| `cycle` | Cycle extracted, e.g. `10_th_cycle` = the stabilized 10th cycle. |

**Example:** `30_30_CPFEM_upto_1_2_percent_10_th_cycle.csv` = the 10th (stabilized) cycle of the CPFEM balanced 30/30 hold simulation at 1.2% strain amplitude.

Inside each CSV, the relevant columns are `strain_zz` (column 33) and `stress_zz` (column 36).

---

## Postprocessing sub-folders

Each sub-folder of *Matlab Code and Data for postprocessing of the results* contains the MATLAB script(s) and CSV data for one analysis:

| Sub-folder | Contents |
|---|---|
| `Accumulated Plastic strain` | Accumulated inelastic strain vs. number of cycles for different hold times. |
| `Damage Map` | Data and script for the damage-partition map (fatigue / mixed / creep regions). |
| `Damage data` | Per-cycle creep and fatigue damage values. |
| `Entropy_Generation and Stress_time plots` | Stress, strain and entropy-generation-rate vs. time over a stabilized cycle. |
| `Fatigue_Test_1%_Validation_Pin_lu` | Hysteresis-loop validation at 1.0% strain (pure fatigue). |
| `Fatigue_Test_1.5%_for_10_cycle_Validation` | Hysteresis-loop validation at 1.5% strain over 10 cycles. |
| `Fatigue_Test_for_1.2_for_10cyle_30_sec_dwell_in_tension` | Creep-fatigue hysteresis at 1.2% strain, 30 s tensile dwell. |
| `Four Hysteris loop` | The four hold-configuration hysteresis loops together. |
| `Hold time comparision` | Effect of tensile hold duration (0–300 s) on response, relaxation and life. |
| `Life at 760 and 980` | Predicted life vs. strain range at 760 °C and 980 °C. |
| `Life data` | Predicted-vs-experimental life (scatter-band comparison). |
| `Pure fatigue and different R ratio effect` | Pure-fatigue response and effect of strain ratio R. |
| `R ratio effect on Life` | Effect of strain ratio R on predicted life. |
| `Temperature_Dependent_Damage_Analysis` | Total damage vs. strain amplitude at 760 °C and 980 °C. |
| `stress_amplitude` | Stress amplitude evolution over cycles (cyclic hardening) vs. experiment. |
| `tension_dwell_1.2%` | Creep-fatigue response at 1.2% strain with a tensile dwell. |

> Rename / re-map these to the exact figure numbers in the final paper as needed.

---

## Damage calculation

The entropy generation rate is computed inside the material subroutine `RVE_Paper_SXNickel.C` and written to the CSV output. The MATLAB postprocessing scripts (in the *Damage data* / *Damage Map* sub-folders) then read this entropy data and evaluate the creep and fatigue damage per cycle, which are combined to give the damage-partition map and the predicted creep-fatigue life.

The governing equations — the entropy generation rate, the fatigue damage, and the creep damage — together with all symbol definitions and material constants, are given in the associated paper (Eqs. 13–15 and Table 1). Refer to the paper for the full formulation.
