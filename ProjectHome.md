## Introduction ##

MOLDY is a parallelised OpenMP short-ranged molecular dynamics program, first written at Harwell Laboratory in the 1980s.  The program is rewritten in a modular fashion to allow for easy user modification, in particular the implementation of new interatomic potentials.  Using Link Cells and Neighbour Lists, the code fully exploits the short range of the potentials, and the slow diffusion expected for solid systems.

The code allows for a wide variety of boundary conditions, including constant pressure, temperature and strain rate.  It also incorporates molecular statics via the conjugate gradients minimisation of the enthalpy.

The code will enable simulation of millions of atoms using short range potentials.  Currently modules for Embedded Atom, Finnis-Sinclair, Lennard Jones and Morse potentials exist.  In addition, the "magnetic" potential formalism of Ackland and Wallenius is available for separate compilation. Alloys containing a number of elements can be simulated, subject only to the available potentials.

## Getting started ##

To get started with using MOLDY head over to the [wiki](https://www.wiki.ed.ac.uk/display/ComputerSim/MOLDY) that contains information on how to compile and run MOLDY as well as some examples.