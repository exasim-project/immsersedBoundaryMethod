# immersedBoundaryMethod

> [!IMPORTANT]
> immsersedBoundaryMethod is under active development!
> If you're interested in contributing to the project, contact the developers or open an [issue](https://github.com/exasim-project/immsersedBoundaryMethod/issues) or a [PR](https://github.com/exasim-project/immsersedBoundaryMethod/pulls)!


## Overview
An OpenFOAM implementation of the Immersed Boundary Method (IBM) developed by

P. Luchini, D. Gatti, A. Chiarini, F. Gattere, M. Atzori & M. Quadrio (2025). _A simple and efficient second-order immersed-boundary method for the incompressible Navier–Stokes equations_, J. Comp. Phys. 539:114245, doi:[10.1016/j.jcp.2025.114245](https://doi.org/10.1016/j.jcp.2025.114245)

The method allows to reproduce a solid boundary within a fluid mesh without the need for meshing: only an STL file representing the solid boundary is required. 

## Compilation

Source OpenFOAM 
```
foamInit
```
Clone the repository
```
git@github.com:exasim-project/immsersedBoundaryMethod.git
```
Enter the repository and initialise the required submodules with
```
cd immsersedBoundaryMethod
git submodules init
git submodules update
```
Source the environment variables required by the repository
```
source bashrc
```
Compile the sources with OpenFOAM's ```wmake```
```
./Allwmake
```

## OpenFOAM Compatibility
This library is developed against and primarily tested with [OpenFOAM v2106](https://www.openfoam.com/news/main-news/openfoam-v2106). 


## License

This project is licensed under the [version 3 of the GNU General Public License](https://www.gnu.org/licenses/gpl-3.0).

For the complete legal text of the license, please refer to the [LICENSE](LICENSE) file located in the root directory of this repository. This file contains the full text of the GNU General Public License.

