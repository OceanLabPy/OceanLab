# OceanLab

Package of Python scripts for Oceanography  (Python +3.6)

## Code Example

Check `examples` folder in our [github repository](github.com/iuryt/OceanLab/examples).

## Installation

`pip install OceanLab`

## Modules

- **OA**
  - *vectoa()*: Objective analysis for vectorial fields;
  - *scaloa()*: Objective analysis for scalar fields;
- **DYN**
  - *dyn_amp()*: Makes the projection of every dynamical mode to velocity to obtain its amplitude;
  - *zeta()*: Calculates the vorticity field by velocity field;
  - *psi2uv()*: Calculates the velocity field by stream function scalar field;
  - *vmodes()*: Calculates the QG pressure modes from N2 profile;
  - *eqmodes()*: Calculates the equatorial pressure and vertical velocity modes from N2 profile;
- **EOF**
  - *eoft()*: Calculates the Empirical Orthogonal Functions;
  - *my_eof_interp()*: Fillgaps on matrix based on EOFs (translated from Cesar Rocha Matlab version);
- **UTILS**
  - *argdistnear()*: Searchs the position of the closest points in an array to a reference point;
  - *meaneddy()*: Performs an eddy-mean decomposing with a low-pass filter;


## Contributors

Everyone can contribute to this code. Some of functions were based on Filipe Fernandes or Cesar Rocha functions and some of them were created with help of Dante C. Napolitano, Hélio M. R. Almeida and Wandrey Watanabe at Ocean Dynamics Lab of University of São Paulo (USP).
