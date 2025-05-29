# PhASER_PDEs

This repository contains the source code associated with the manuscript

Marchi, Khalek, George, Weitz, Chait: Stable coexistence and transport of lytic phage infections with migrating bacterial hosts

It provides the code to run the model and generate the numerical data presented in the paper, and it allows reproduction of all figures of the manuscript from such pre-computed numerical data. 

## Dependencies installation and repository organization 

### Installation requirements

The code uses Python 3.
A number of standard scientific python packages are needed for the numerical simulations and visualizations. The main packages needed to run the code are:
- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/stable/index.html)

An easy way to manage all of the code dependencies is to install a Python distribution such as [Anaconda](https://www.anaconda.com/). All required packages are available in Anaconda. Installing all of them should take less than an hour on a standard desktop machine.

### Files organization/running the code
In [the bash scripts folder](./bash_scripts) there are bash handlers used to sweep through the parameters space and generate movies from pre-processed png frames. The [bash wrapper](./bash_scripts/model_batch_runs.sh) sets two of the model parameters, and calls the python [simulation code](./scripts_generate_figures) running the PDE integration. The parameters needed to generate all manuscript figures are thoroughly listed in the manuscript's Supplementary Text, and are also reported in the simulation code as commented python dictionaries. Note that the code can take 5-20 hours to run depending on the chosen parameters, and the output files typically weight 2-10 GB. 

[The plots folder](./plots) contains the python scripts to reproduce the figures panels presenting numerical results as they are presented in the manuscript, once the corresponding synthetic data are generated. Each script name indicates which figure panel it corresponds to. Figure panels can then be assembled through an image editor such as inkscape. All scripts require specific relative folder placements to retrieve the correct synthetic data. The synthetic data plotted in the manuscript are too heavy to be uploaded to this repository, but all code and information necessary to re-generate the model outputs are provided here and in the manuscript's Supplementary Text.

[Lib](./lib) contains some plotting cosmetics definitions, and the [PDE integration functions](./lib) which implements the Finite Differences algorithm detailed in the manuscript's Supplementary Text.

Finally, [the data folder](./growth_data) contains four csv files with the intensity profile data used to infer the bacteria parameters for all three _E. coli_ strains presented in the manuscript. Specifically, [two independent profiles](./growth_data/Intensity_F8_24cm) were used to infer the F8 strain parameters, and [other two measured profiles](./growth_data/Intensity_F9F10) informed F9 and F10 strains parameters. The two jupyter notebooks contain the code used to process these data and infer the bacteria growth parameters for [F8](./growth_data/import_data_radialCFU_F8fin_clean.ipynb) and [F9/F10](./growth_data/import_data_radialCFU_F9_clean.ipynb) strains respectively. These notebooks generate panels composing Supplementary Fig S11 and S12. 
The [Microsoft Excel file](./growth_data/Scanner_OD.xlsx) contains the calibration data used throughout the analysis to transform OD to CFUs.
