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
| `Data and MATLAB Scripts for Postprocessing` | MATLAB postprocessing scripts and the CSV data used to generate the figures. |
| `README.md` | This file. |

**Source and header files**

| File | Description |
|---|---|
| `RVE_Paper_SXNickel.C` / `.h` | The custom MOOSE material (UMAT-equivalent stress update). Implements the flow rule on the 12 FCC slip systems, dislocation-density evolution, Armstrong–Frederick back-stress with dynamic and Arrhenius static recovery, the implicit local Newton–Raphson update, and the algorithmic tangent. Registered in the app as object type `RVE_Paper_SXNickel`. |

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



# Data and MATLAB Scripts for Postprocessing

This repository contains the MATLAB post-processing scripts and the CSV data
files used to generate the figures in the paper. Each sub-folder corresponds
to one (or a small group of) figures. The source data are outputs of a
dislocation density-based crystal plasticity (CPFEM) model implemented in the
MOOSE framework.

## How to use

1. Open the sub-folder for the figure you want to reproduce.
2. Open the `.m` script inside it in MATLAB.
3. Run the script — it reads the accompanying `.csv` files in the same folder
   and produces the figure.

Every script header names the paper figure it generates
(e.g. `%%%%%% Figure-5(c) %%%%%%`).

## Data file naming convention

The CSV data files use a compact code for hold configuration and strain range.
The most common pattern is:

```
<hold_config>_<strain_range>.csv
```

| Field          | Meaning |
|----------------|---------|
| `hold_config`  | Dwell configuration as `<tension>_<compression>` hold seconds: `0_0` (pure fatigue, no hold), `30_0` / `60_0` / `120_0` / `240_0` / `300_0` (tensile hold), `0_60` (compressive hold), `30_30` (balanced tension–compression hold). |
| `strain_range` | Total strain range Δε in %, with an underscore for the decimal point (e.g. `2_4` = 2.4 % range → strain amplitude εₐ = 1.2 %; `1_6` = 1.6 % → εₐ = 0.8 %). |

**Strain-code → amplitude reference** (amplitude εₐ = ½ × range):

| Code | Range Δε | Amplitude εₐ |
|------|----------|--------------|
| `1_4` | 1.4 % | 0.7 % |
| `1_6` | 1.6 % | 0.8 % |
| `1_8` | 1.8 % | 0.9 % |
| `2_0` | 2.0 % | 1.0 % |
| `2_2` | 2.2 % | 1.1 % |
| `2_4` | 2.4 % | 1.2 % |

**Example:** `30_30_2_4.csv` = balanced 30/30 s hold simulation at 2.4 % total
strain range (1.2 % amplitude).

Some files carry an extra descriptive prefix or suffix instead of the bare code
— for example `Fatigue_60_0_2_4.csv` (fatigue / creep-fatigue set),
`CFI_2_4_30_30_new.csv` (creep-fatigue interaction, strain range then hold),
and the `shi_…`, `hold_…`, `cyclic_…`, `Pin_…` variants. Files beginning with
`R_…` are strain-ratio R cases, and `mean_…` / `mean_strain_…` are mean-strain
cases. A trailing `_new` or `_for_curve` marks a revised or curve-fitting
version of the same case.

Inside each CSV, the two most relevant columns are `strain_zz` (**column 33**)
and `stress_zz` (**column 36**).

## Folder-to-figure map

| Figure        | Folder                                          | Script                                   | Description |
|---------------|-------------------------------------------------|------------------------------------------|-------------|
| Fig. 2(a)     | `Cyclic_test_Hysteresis_loop`                   | `Both_1_1_5_percentage.m`                | Pure-fatigue hysteresis at 1.0 % and 1.5 % strain |
| Fig. 2(b)     | `Cyclic_test_Hysteresis_loop`                   | `Tension_test_dwell_60sec_in_tension.m`  | 60 s tensile-dwell hysteresis loop |
| Fig. 2(c)     | `Cyclic_test_Hysteresis_loop`                   | `Cyclic_60_Hold_Compression.m`           | 60 s compressive-hold hysteresis loop |
| Fig. 2(d)     | `Cyclic_test_Hysteresis_loop`                   | `Creep_fatigue_test.m`                   | Combined creep-fatigue hysteresis loop |
| Fig. 3        | `Four Hysteresis loop`                          | `four_loop_hysteresis.m`                 | The four hold-configuration hysteresis loops together |
| Fig. 4(a)     | `Life at 760 and 980/Figure_4_(a)_paper`        | `combined_plot.m`                        | Predicted life vs. strain range at 760 and 980 °C |
| Fig. 4(b)     | `stress_amplitude`                              | `stress_amp.m`                           | Stress-amplitude evolution over cycles (cyclic hardening) vs. experiment |
| Fig. 5(a,b)   | `Entropy_Generation and Stress_time plots`      | `stress_strain_time.m`                   | Stress/strain and entropy-generation rate vs. time over a stabilized cycle |
| Fig. 5(c)     | `Damage creep and fatigue evolution`            | `Figure_5_c.m`                           | Creep and fatigue damage vs. time (30/30, 1.2 %) |
| Fig. 6(a)     | `Predicted_vs_Experimental_life`                | `Figure_6_a.m`                           | Predicted-vs-experimental life (scatter-band comparison) |
| Fig. 6(b)     | `Predicted_vs_Experimental_life`                | `Figure_6_b.m`                           | Predicted-vs-experimental life (scatter-band comparison) |
| Fig. 7(a–d)   | `Hold time comparison`                          | `three_time_60_120_240s.m`               | Effect of tensile hold duration (60/120/240 s) at 980 °C |
| Fig. 8        | `Accumulated Plastic strain`                    | `all_plastic_accumulation.m`             | Accumulated inelastic strain vs. number of cycles |
| Fig. 9(a)     | `Damage data`                                   | `corrected_damage.m`                     | Per-cycle creep and fatigue damage values |
| Fig. 10       | `Hold time comparison`                          | `map_2_new.m`                            | Hold-time comparison / damage-partition map |
| Fig. 11(a,b)  | `Mean strain plot`                              | `Figure_11_a_b.m`                        | Mean-strain plots |
| Fig. 12(a)    | `Pure fatigue and different R ratio effect`     | `No_Hold_Fatigue.m`                      | Pure-fatigue response and effect of strain ratio R |
| Fig. 12(b)    | `R ratio effect on Life`                        | `Damage_comparison.m`                    | Effect of strain ratio R on predicted life |
| Fig. 12(c)    | `strain_ratio_life`                             | `Figure_12_c.m`                          | Effect of strain ratio R on predicted life (strain-ratio view) |
| Fig. 13(a)    | `Life data`                                     | `Figure_13_a.m`                          | Life data, 30/30 configuration |
| Fig. 13(b)    | `Life data`                                     | `Figure_13_b.m`                          | Life data, 60/0 configuration |
| Fig. 14(a)    | `Life at 760 and 980`                           | `Figure_14_a.m`                          | Predicted life vs. strain range at 760 °C |
| Fig. 14(b)    | `Life at 760 and 980`                           | `Figure_14_b.m`                          | Predicted life vs. strain range at 980 °C |
| Fig. 15       | `Temperature_Dependent_Damage_Analysis`         | `new_temp_depen.m`                       | Total damage vs. strain amplitude at 760 and 980 °C |

**Supporting folder (not a numbered paper figure):** `Damage Map` — grid data
and script used to build the fatigue / mixed / creep damage-partition map.

## Requirements

- MATLAB (tested on R2025b).
- No additional toolboxes are required for the basic plots.

## Citation

If you use these scripts or data, please cite the associated paper.
---

## Damage calculation

The entropy generation rate is computed inside the material subroutine `RVE_Paper_SXNickel.C` and written to the CSV output. The MATLAB postprocessing scripts (in the *Damage data* / *Damage Map* sub-folders) then read this entropy data and evaluate the creep and fatigue damage per cycle, which are combined to give the damage-partition map and the predicted creep-fatigue life.

The governing equations — the entropy generation rate, the fatigue damage, and the creep damage — together with all symbol definitions and material constants, are given in the associated paper (Eqs. 13–15 and Table 1). Refer to the paper for the full formulation.
