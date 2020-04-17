Our approach is based on the interface method and aims to calculate the melting temperature of materials fully automatically.

The interface method is a well established approach for predicting melting points of materials with interatomic potentials. 
However, applying the interface method is tedious and involves significant human intervention, even in the case when fast 
interatomic potentials are employed. The whole procedure involves several successive tasks: estimate a rough melting point, 
set up the interface structure, run molecular dynamic calculations and analyse the data. Loop calculations are necessary if 
the predicted melting point is different from the estimated one by more than a certain convergence criterion, or if full 
melting/solidification occurs. In this case monitoring the solid-liquid phase transition in the interface structure becomes 
critical. As different initial random seeds for the molecular dynamic simulations within the interface method induce slightly 
different melting points, a few ten or hundred interface method calculations with different random seeds are necessary for 
performing a statistical analysis on these melting points. Considering all these technical details, the work load for 
manually executing and combining the various involved scripts and programs quickly becomes prohibitive. To simplify and 
automatize the whole procedure, we have implemented the interface method into pyiron (http://pyiron.org).

The fully automatized procedure allows to efficiently and precisely predict melting points of unaries represented by arbitrary potentials with only two user-specified parameters (potential and element).

To perform the automatized interface method calculations, the Jupyter notebook needs to be downloaded. The notebook can be executed after modifying the potential and element. Running directly the notebook the details of each step are automatically plotted and can be analysed retrospectively. To further simplify the application of our tool, we provide another approach. This requires an input file, which contains the two user-specified parameters, element and interatomic potential, in the following format:

{
  "config": [ "pair_style eam/alloy \n",
              "pair_coeff * * potential Al\n"
             ],
  "filename": ["./example/potential"],
  "species": ["Al"],
  "element": "Al"
}

Here ”config” provides the LAMMPS parameters for loading the potential, ”filename” shows the location of the potential file, ”species” gives the species implemented in the potential, and ”element” defines on which material the melting point calcula- tions will be performed. The input file is in a JSON (JavaScript Object Notation) format that is easy to read and write for hu- mans and easy to parse and generate for machines. After creat- ing the input.json file, one can simply execute the protocol by using snakemake:

snakemake --use-conda

The parameters defined in the input.json file will overwrite those in the Jupyter notebook. Using this approach, there is no need to interfere with all the computational and technical details. Both approaches write the estimated melting points and predicted melting points from all loop calculations in a file melting.json.
