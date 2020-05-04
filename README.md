# Introduction
The [pyiron](http://pyiron.org) based module and Jupyter Notebook *Melting* allow the fully automated computation of melting points of unary crystals for arbitrary interatomic potentials that are compatible with the molecular dynamics engine [LAMMPS](https://lammps.sandia.gov). It is based on the interface method where the evolution of the solid and the liquid phase are monitored as function of temperature. The only mandatory input parameters required are the chemical element and the interatomic potential file. 

# Use the melting point method in pyiron
[pyiron](http://pyiron.org) is an integrated development environment (IDE) for computational materials science. It was used to develop the melting point protocol and is designed for the development of complex simulation protocols in general. If you are already a pyiron user or want to understand the detailed steps of the melting point method we recommend using the melting point protocol with pyiron. If you are only interested in the calculated melting point values, the snakemake approach explained below might be more suitable for you. Both approaches are limited to the linux operation system and have been successfully tested with the linux subsystem for windows. 

## Installation 
Please install the following packages: 

- [pyiron](http://pyiron.org)
- [Ovito](https://ovito.org)
- [Lammps](https://lammps.sandia.gov)

All packages are available via conda-forge and can be installed with the following command: 
```
conda install -c conda-forge pyiron nglview lammps jupyter scikit-learn ovito mscorefonts
```

For the installation of pyiron and the configuration of Lammps and NGLview within pyiron please refer to the [pyiron manual](https://pyiron.readthedocs.io/en/latest/source/installation.html).

## Run the Jupyter Notebook
You can directly download the Jupyter notebook from this Github repositorry [*script.ipynb*](https://github.com/pyrion/pyrion_meltingpoint/blob/master/scripts/script.ipynb) copy it to your pyiron projects folder and execute it there. In line 10 and 11 the input parameters can be modified to select a custom potential and change parameters of the melting point calculation. After the calculation finished successfully it creates an *output.json* file which contains the final melting point prediction as well as the intermediate results. 

## Analyse 
You can analyse the melting point calculation using the [*plot.ipynb*](https://github.com/pyrion/pyrion_meltingpoint/blob/master/scripts/plot.ipynb) notebook. 

# Use the melting point method with snakemake 
In contrast to the pyiron approach snakemake drastically reduces the number of input parameters available to the user. Snakemake handles the setup of pyiron and the installation of the dependencies via condda. This method is recommended for high throughput calculation as well as automated validation of interatomic potentials. Still it might also be sufficient for users who just want to calculate the melting point for a given interatomic potential. Both approaches are limited to the linux operation system and have been successfully tested with the linux subsystem for windows. 

Start with installing snakemake from conda: 
```
conda install -c bioconda -c conda-forge snakemake=5.16
```

Copy one of the *input.json* files from the examples to the root of this directory: 
```
cp examples/bccFe/input.json .
```

The content of the *input.json* is: 
```
{
    "config": [
        "pair_style eam/alloy \n", 
        "pair_coeff * * Fe-C-Bec07.eam Fe C\n"
    ], 
    "filename": "./examples/bccFe/Fe-C-Bec07.eam", 
    "species": ["Fe", "C"], 
    "element": "Fe"
}
```
After copying the *input.json* it can be executed using: 
```
snakemake --use-conda --cores 1 
```
The parameters defined in the *input.json* file will overwrite those in the Jupyter notebook. With this approach, there is no need to interfere with all the computational and technical details. 
    
The results are saved in the *output.json* file and can be analysed with the [*plot.ipynb*](https://github.com/pyrion/pyrion_meltingpoint/blob/master/scripts/plot.ipynb) notebook. 

# FAQ
## How to run in parallel? 
A single melting point calculation takes 50-100 CPU hours, so it makes a lot of sense to run the code in parallel. While the protocol itself is written in a serial way, the individual Lammps calculation can be executed in parallel. To enable parallel execution inset the option: 
```
"cpu_cores": 8,
```
Either in line 10 of the jupyter notebook or in the *input.json* file. When snakemake is used it is not necessary to increase the `--cores` count in the snakemake command. 

# Acknowledgments
If you use the melting point protocol in your scientific work, please consider citing:
```
  @article{melting-paper,
    title = {A fully automated approach to determine the melting temperature of crystalline materials},
    author = {Li-Fang Zhu and Jan Janssen and Shoji Ishibashi and Blazej Grabowski and Jörg Neugebauer},
    journal = {in preparation}
  }

  @article{pyiron-paper,
    title = {pyiron: An integrated development environment for computational materials science},
    journal = {Computational Materials Science},
    volume = {163},
    pages = {24 - 36},
    year = {2019},
    issn = {0927-0256},
    doi = {https://doi.org/10.1016/j.commatsci.2018.07.043},
    url = {http://www.sciencedirect.com/science/article/pii/S0927025618304786},
    author = {Jan Janssen and Sudarsan Surendralal and Yury Lysogorskiy and Mira Todorova and Tilmann Hickel and Ralf Drautz and Jörg Neugebauer},
    keywords = {Modelling workflow, Integrated development environment, Complex simulation protocols},
  }
```
