# Introduction
The [pyiron](http://pyiron.org) based module and Jupyter Notebook *Melting* allow the fully 
automated computation of melting points
of unary crystals for arbitrary interatomic potentials that are compatible with the molecular dynamics engine 
[LAMMPS](https://lammps.sandia.gov). It is based on the interface method where the evolution of the solid and 
the liquid phase are monitored as function of temperature. The only mandatory input parameters required 
are the chemical element and the interatomic potential file. 

# Set up the environment
To run the module *Melting*, a linux system is recommended. Under Windows please install the Linux subsystem. 
Please install:

- [pyiron](https://pyiron.github.io/source/installation.html)
- [Ovito](https://anaconda.org/conda-forge/ovito) (To run *Melting* in a python console)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) (Optional: To 
provide input to the notebook without having to edit it.)

   
# Run the Jupyter Notebook
Download the Jupyter notebook, [*melting.ipynb*](https://github.com/lfzhu-phys/Melting/blob/master/melting.ipynb), 
to the path where the pyiron projects are located, 
e.g., "~/pyiron/projects/". Two approaches can be used to run it:

- Approach 1:
    ```python
    cd ~/pyiron/projects
    jupyter notebook
    ```

    Open the downloaded notebook, modify the potential and element, and execute the notebook directly.

- Approach 2:

    Use snakemake to run the notebook. This approach requires an *input.json* file in the following format:

    ```json
    {
     "config":   [ "pair_style eam/alloy \n",
                   "pair_coeff * * potential Al\n"],
     "filename": ["/home/pyiron/projects/potential"],
     "species":  ["Al"],
     "element":  "Al",
    }
    ```
    and a Snakefile, which is in the following format:
    ```json
    rule melting:
        notebook:
            "melting.ipynb"
    ```
    After creating *input.json* and *Snakefile* in the same directory as the notebook *melting.ipynb*, 
    execute the protocol by:
    ```
    snakemake --use-conda 
    ```
    The parameters defined in the input.json file will overwrite those in the Jupyter notebook. With this approach, there is no need to interfere with all the computational and technical details. 
    
Both approaches write the estimated melting points and predicted melting points from all loop calculations 
in a file *melting.json*.

# Acknowledgments
If you use the melting point protocol in your scientific work, please consider citing:
```
  @article{melting-paper,
    title = {A fully automated approach to determine the melting temperature of crystalline materials},
    author = {Li-Fang Zhu and Jan Janssen and Shoji Ishibashia and Blazej Grabowskic and Jörg Neugebauer},
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
