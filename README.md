# PhASER_PDEs

This repository contains the source code associated with the manuscript

Marchi, Khalek, George, Weitz, Chait: Stable coexistence and transport of lytic phage infections with migrating bacterial hosts

It allows reproduction of all figures of the manuscript.

## Figures reproduction code

#### Installation requirements

The visualization code uses Python 3.

A number of standard scientific python packages are needed for the numerical simulations and visualizations. An easy way to install all of these is to install a Python distribution such as [Anaconda](https://www.anaconda.com/):

- [numpy](https://numpy.org/)
- [scipy](https://www.scipy.org/)
- [matplotlib](https://matplotlib.org/stable/index.html)


#### Files organization/running the code
In  [the bash scripts folder](./scripts_generate_figures) there are the bash handlers used to sweep through parameters space and call the ODE models integrators (in [models](./models_python)), generating synthetic data and organizing the folders in the required relative positions. The scripts names are self-explanatory. [The plots folder](./plots) contains the scripts to produce the figures once synthetic data are generated.
[Lib](./lib) contains some plotting cosmetics definitions. Finally [the jupyter notebook](./data_processing_scripts/jupyter-notebook/import_data_LD_2023_clean.ipynb) contains the data processing code that takes Microsoft Excel files (in [the data folder](./data) ) with the fluctuation test results and produces txt files compatible with the inference tool [bzrates](http://www.lcqb.upmc.fr/bzrates#:~:text=bz%2Drates%20is%20a%20web,mutation%20rates%20from%20fluctuation%20assays.).
