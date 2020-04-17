Our approach is based on the interface method and aims to calculate the melting temperature of materials fully automatically.

The interface method is a well established approach for predicting melting points of materials with interatomic potentials. 
How- ever, applying the interface method is tedious and involves significant human intervention, even in the case when fast 
interatomic potentials are employed. The whole procedure involves several successive tasks: estimate a rough melting point, 
set up the interface structure, run molecular dynamic calculations and analyse the data. Loop calculations are necessary if 
the predicted melting point is different from the estimated one by more than a certain convergence criterion, or if full 
melting/solidification occurs. In this case monitoring the solid-liquid phase transition in the interface structure becomes 
critical. As different initial random seeds for the molecular dynamic simulations within the interface method induce slightly 
different melting points, a few ten or hundred interface method calculations with different random seeds are necessary for 
performing a statistical analysis on these melting points. Considering all these technical details, the work load for 
manually executing and combining the various involved scripts and programs quickly becomes prohibitive. To simplify and 
automatize the whole procedure, we have implemented the interface method into pyiron (http://pyiron.org).
