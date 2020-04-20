# Introduction
Melting is based on the interface method and aims to fully automatically calculate the melting temperature of unaries represented by arbitrary interatomic potentials. With this approach, only two parameters, element and potential, are necessary to specify. The automated procedure is fully implemented in a pyiron-based Jupyter notebook (http://pyiron.org) and LAMMPS is the engine to perform all the automatic calculations.

# Set up the environment
To perform Melting, linux system is recommended. The user needs to install:


1) pyiron, NGLview visualization framework, and LAMMPS.


   Guide for installation, see: https://pyiron.github.io/source/installation.html

2) Ovito.


   Guide for installation, see: https://anaconda.org/conda-forge/ovito
   
   
in case you want to run the notebook in a terminal:


3) snakemake.


   Guide for installation, see: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html
# Run the notebook
The user needs to download the Jupyter notebook, melting.ipynb, to the path where the pyiron projects locate, e.g., "~/pyiron/projects/". Two approaches can be used for running the notebook.

Approach 1:
```python
cd ~/pyiron/projects
jupyter notebook
```

Open the downloaded notebook, modify the potential and element, and execute the notebook directly.

Approach 2:

Use snakemake to run the notebook. This approach requires an input.json file, which contains the two user-specified parameters, element and interatomic potential, in the following format:

```json
{
 "config":   [ "pair_style eam/alloy \n",
               "pair_coeff * * potential Al\n"],
 "filename": ["~/pyiron/projects/potential"],
 "species":  ["Al"],
 "element":  "Al"
}
```
and a Snakefile, which is in the following format:
```json
rule melting:
    input:
        "input.json"
    output:
        "melting/melting.json
    notebook:
        "melting.ipynb"
```
After creating input.json and Snakefile, one can simply execute the protocol by:
```python
snakemake --use-conda --snakefile Snakefile --cores n
```
The parameters defined in the input.json file will overwrite those in the Jupyter notebook. With this approach, there is no need to interfere with all the computational and technical details. 

Both approaches write the estimated melting points and predicted melting points from all loop calculations in a file melting.json.
