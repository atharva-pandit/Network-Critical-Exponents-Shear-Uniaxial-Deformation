## Network-Critical-Exponents-Shear-Uniaxial-Deformation
Simulations of 2D disordered fibre networks. Apply shear and uniaxial compression or deformation. Measure critical exponents.

### Dependencies

* C++20 or higher
* Python 3.10 or higher

### Execution

* network.cpp - Generate network of size W, connectivity Z & version V.
* shear.cpp - Pre-stress the network to deformation E. Apply incremental shear.
* shear_spring.cpp - Same as shear.cpp, but without bending interactions.
* nonaffinity.cpp - Measure non-affine rearrangements from sheared and relaxed states.
* setup.hpp - Define structs and basic functions.
* read_write.hpp - Read/write intermediate states or LAMMPS compatible input/dump files for visualisation.
* minimise.hpp - Energy minimisation & MD integration.
* minimise_spring.hpp - Same as minimise.hpp, but without bending interactions.
* slide.ipynb - Plot data & determine critical exponents.

* Run programs with optimisations turned on & provide network information as input from terminal.
```
g++ shear.cpp -O3 -march=native
./a.out 100 3.2 6 -0.1 0 1
```

## Author
Atharva Pandit - atharva.pandit@uni-a.de

## Version History
* 0.1
    * Initial Release

## License

This project is licensed under the GNU General Public License - see the LICENSE.md file for details

## Acknowledgments
I would like to acknowledge you for bothering to read this far.
