## README
This repository contains codes for simulations of 2D disordered fibre networks, presented in the paper "**Strain-stiffening critical exponents of fiber networks under uniaxial deformation**".

## Dependencies
* C++20 or higher
* Python 3.10 or higher

## Codes
* **network.cpp** - Generate a network of size W, connectivity Z & version V.
* **shear.cpp** - Pre-stress the network to deformation E. Apply incremental shear in direction R.
* **shear_spring.cpp** - Same as shear.cpp, but without bending interactions.
* **nonaffinity.cpp** - Measure non-affine rearrangements from sheared and relaxed states.
* **setup.hpp** - Define structs and basic functions.
* **read_write.hpp** - Read/write intermediate states or LAMMPS compatible input/dump files for visualisation.
* **minimise.hpp** - Energy minimisation, MD integration & stress calculation.
* **minimise_spring.hpp** - Same as minimise.hpp, but without bending interactions.
* **slide.ipynb** - Plot data & determine critical exponents.

## Author
**Atharva Pandit**\
Contact: atharva.pandit@uni-a.de

## Version History
* 0.1
    * Initial Release

## License
This project is licensed under the GNU General Public License v3.0 - see the LICENSE.md file for details

## Acknowledgments
I would like to acknowledge you for bothering to read this far.
