Melting is based on the interface method and aims to fully automatically predict the melting temperature of unaries represented by arbitrary interatomic potentials. With this approach, only two parameters, element and potential, are necessary to specify. The automated procedure is realized by implementing the full simulation protocal into a pyiron-based Jupyter notebook. (http://pyiron.org).

To perform Melting (linux system is recommended), the user needs:

1) install pyiron, NGLview visualization framework, and Lammps.
   Guide for installation see: https://pyiron.github.io/source/installation.html

2) install ovito.
   Guide for installation see: https://anaconda.org/conda-forge/ovito

In case you want to run the notebook in a terminal:

3) install snakemake.
   Guide for installation see: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

After setting up the environment, the user needs to download the Jupyter notebook to the path of the pyiron project, e.g., "~/pyiron/projects/", and go to the path. Then two approaches can be used for running the notebook.

Approach 1:

cd ~/pyiron/projects
jupyter notebook

Open the downloaded notebook, modify the potential and element, and execute the notebook directly.

Approach 2:

Use snakemake to run the notebook. This approach requires an input.json file, which contains the two user-specified parameters, element and interatomic potential, in the following format:

{
 "config":   [ "pair_style eam/alloy \n",
               "pair_coeff * * potential Al\n"],
 "filename": ["./example/potential"],
 "species":  ["Al"],
 "element":  "Al"
}

After creating the input.json file, one can simply execute the protocol by using snakemake:

snakemake --use-conda

The parameters defined in the input.json file will overwrite those in the Jupyter notebook. 

With both approaches, there is no need to interfere with all the computational and technical details. Theestimated melting points and predicted melting points from all loop calculations will be written in a file melting.json.
