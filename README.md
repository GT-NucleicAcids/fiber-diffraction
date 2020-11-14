[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/GT-NucleicAcids/fiber-diffraction/HEAD?filepath=diffraction.ipynb)

# Overview
The repository contains a script for simulating fiber diffraction pattern for molecular systems. It takes the 3D geometry of a molecular system and generates the fiber diffraction pattern.

The script can generate model helical noncovalent polymers given the 3D geometry of the monomeric unit, the rise, the twist, and the number of units in the stack. The repository contains example geometries for hexameric rosette units based on triaminopyrimidine and cyanuric acid. The examples can be run using Jupyter notebook. Click on the `launch binder` icon to try the script online.

# Files
* `diffraction.ipynb`: A Jupyter notebook containing example runtime options.
* `scripts.py`: A python file containing the functions used for simulating the diffraction pattern.
* `geometries/TAP_Cy.xyz`: The 3D geometry of a hexameric rosette unit.
* `geometries/TAP_4MCyCo6.xyz`: The 3D geometry of a hexameric rosette unit with chiral exocyclic tails.
* `geometries/TAP_CyCo6.xyz`: The 3D geometry of a hexameric rosette unit with achiral exocyclic tails.