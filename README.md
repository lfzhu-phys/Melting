# Introduction
Melting is based on the interface method and aims to fully automatically predict the melting temperature of unaries represented by arbitrary interatomic potentials. With this approach, only two parameters, element and potential, are necessary to specify. The automated procedure is fully implemented in a [pyiron](http://pyiron.org)-based Jupyter notebook and [LAMMPS](https://lammps.sandia.gov) is the engine to perform all the automatic calculations.

# Set up the environment
To perform Melting, linux system is recommended. The user needs to install:


1) pyiron, NGLview visualization framework, and LAMMPS.


   Guide for installation, see: https://pyiron.github.io/source/installation.html

2) Ovito.


   Guide for installation, see: https://anaconda.org/conda-forge/ovito
   
   
in case you want to run the notebook in a terminal,


3) snakemake.


   Guide for installation, see: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
# Run the notebook
The user needs to download the Jupyter notebook, melting.ipynb, to the path where the pyiron projects locate, e.g., "~/pyiron/projects/". Two approaches can be used to run the notebook.

Approach 1:
```python
cd ~/pyiron/projects
jupyter notebook
```

Open the downloaded notebook, modify the potential and element, and execute the notebook directly.

Approach 2:

Use snakemake to run the notebook. This approach requires an input.json file in the following format:

```json
{
 "Config":   [ "pair_style eam/alloy \n",
               "pair_coeff * * potential Al\n"],
 "Filename": ["/home/pyiron/projects/potential"],
 "Species":  ["Al"],
 "element":  "Al",
 "crystalstructure": "Fcc"
}
```
and a Snakefile, which is in the following format:
```json
rule melting:
    notebook:
        "melting.ipynb"
```
After creating input.json and Snakefile in the same path as that melting.ipynb locates, one can simply execute the protocol by:
```
snakemake --use-conda 
```
The parameters defined in the input.json file will overwrite those in the Jupyter notebook. With this approach, there is no need to interfere with all the computational and technical details. 

Both approaches write the estimated melting points and predicted melting points from all loop calculations in a file melting.json.
