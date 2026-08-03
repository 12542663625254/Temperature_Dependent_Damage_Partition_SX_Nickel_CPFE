# Temperature_Dependent_Damage_Partition_SX_Nickel_CPFE
Dislocation density-based crystal plasticity (CPFE) model for temperature-dependent creep-fatigue damage partitioning in single-crystal nickel superalloys. MOOSE implementation with Arrhenius static-recovery and an explicit temperature-dependent creep-fatigue interaction term. The code and the data are attached here.



##### HOW TO RUN THE CODE

## Running the simulations in MOOSE

### Prerequisites
- A working MOOSE installation (see https://mooseframework.inl.gov/getting_started/installation/)
- The custom material model in `src/` and `include/` compiled into your MOOSE app
- PETSc (bundled with MOOSE)

### Building the custom model
The custom crystal plasticity material (stress update, local Newton–Raphson,
and algorithmic tangent) is in `src/materials/` with headers in `include/materials/`.
Copy these into your MOOSE app's source tree and rebuild:

    make -j4

### Running an input file
Each case has a MOOSE input file (`.i`) in the corresponding folder. Run with:

    ./your-app-opt -i "Input file for MOOSE/<case>.i and more precisely it is given in Different Loading condition and change the input file according to strain amplitude and hold configurations"

Replace `your-app-opt` with the name of your compiled MOOSE executable.

### Folder guide
- `src/`, `include/`  — custom material model source and headers
- `Input file for MOOSE/`      — MOOSE `.i` input files
- `Different Loading condition/` — input variants for each loading case (0/0, 60/0, 0/60, 30/30) and change the R-ratio effect as mean strain which can be done by input file, the temperature effect for 760°C and 980°C can be calculated manually and do the changes in the fcc_props.in file accordingly)
- `Matlab Code and Data for Post processing of the results/` — post-processing scripts and figure data
