# Q_Poisson1D

This code solves the Poison equation in 1D on inhomogeneous grid in semiconductors heterostructures. As a results, it gives the band bending profile for any heterostructures. It is specially useful for CMOS, diode and whatever npn & pnp junction. 
A strain model is included. It basically shifts the conduction and valence band edge
The code does NOT take into account the non-parabolicity of the bands. The electron and the hole masses are loaded from the "material.csv" file for each material.
Additional material can be added in the "materialDB_ZB.csv" file

![image](https://user-images.githubusercontent.com/35040499/113292239-866e8e80-92f4-11eb-8791-5bf75c9489a8.png)

