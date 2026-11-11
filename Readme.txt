

This repository provides the MATLAB source code for the paper: 
**"Structurally targeted immunization strategy for suppressing epidemic spreading in directed networks"**

## Requirements
- MATLAB (R2020a or later recommended)
- MATLAB Bioinformatics Toolbox (for `graphconncomp` support, or use `conncomp` in newer versions)

## Directory Structure
- `main_real_network.m`: Main script to run the real-world network experiments.
- `src/`: Source code folder.
  - `directed_SI_simulation/`: SI epidemic models.
  - `immunization_strategy/`: Implementation of HDA, CI, TAD, and DTS- (GCC/DHI).
- `data/`:'Real_Network_SI_data.mat'.


