# SVEIHRM_model

These codes are experimental reproductions of the paper "A model of COVID-19 pandemic with vaccines and mutant viruses" published in PLOS one.

1. The SVEIHRM_codes folder contains code and data files (.csv) to reproduce the figures from Figures 3 to 15.
We can get the result by opening and running the desired figure file, but our code includes matlab's paid package (optimization tool).

- The diffun_m.m file represents our ODE system.
- The ParaToODE.m file is a code that estimates parameters according to the ODE system set in diffun_m.m.
- These two files are used in different figure codes.
- Files with the .mat extension are stored variables. These variables can be saved again using figure codes.

2. The resulting figures from the code are attached to the Figures folder.
These figures are saved as tif files.
If you want to save figures, run the code and save the resulting figure with the desired extension.
