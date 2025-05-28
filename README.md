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
In [the bash scripts folder](./scripts_generate_figures) there are bash handlers used to sweep through the parameters space and generate movies from pre-processed png frames. The [bash wrapper](./scripts_generate_figures) sets two of the model parameters, and calls the python [simulation code](./scripts_generate_figures) running the PDE integration. The parameters needed to generate all manuscript figures are thoroughly listed in the manuscript's Supplementary Text, and are also reported in the simulation code as commented python dictionaries. Note that the code can take 5-20 hours to run depending on the chosen parameters, and the output files tipically weight 2-10 GB. 

[The plots folder](./plots) contains the python scripts to produce the figures once synthetic data are generated. They reqire specific relative folder placements to find the correct data. Interested users should feel free to adapt the input data paths in each script to their local files organization. 

[Lib](./lib) contains some plotting cosmetics definitions, and the [PDE integration functions](./lib). 

Finally [the jupyter notebook](./data_processing_scripts/jupyter-notebook/import_data_LD_2023_clean.ipynb) contains the data processing code that takes Microsoft Excel files (in [the data folder](./data) ) with the fluctuation test results and produces txt files compatible with the inference tool [bzrates](http://www.lcqb.upmc.fr/bzrates#:~:text=bz%2Drates%20is%20a%20web,mutation%20rates%20from%20fluctuation%20assays.).
