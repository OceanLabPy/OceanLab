# OceanLab

Package of Python scripts for Oceanography - Version 1.0

## Synopsis

This lib was made by **Iury Sousa** (<mailto:simoesiury@gmail.com>) from Brazil to help him and other people that
needs to organize, process and plot CTD and VADCP data from cruise. Mainly directioned
to mesoscale oceanographers. The number of functions and subpackages and uses increased in process of time, and now
it's not restricted to mesoscale physical oceanographers.

Now this package have this subpackages below:

- **ADCP**
- **AO**
- **CTD**
- **DYN** 
- **EOF**
- **SEAPLOT**
- **utils**


## Code Example

TODO

## Motivation

There's still a lack of packages and tutorials for Python uses on oceanography. So, the OceanLab
project is based not just provide a package, but to create a teaching tool for oceanographers though.
The website with Python instalation and tutorials is still only available on **Portuguese**, but I'm 
working on english version. See the website on this link: [OceanLab][1]. Help me to 

## Installation

Clone this directory and run:

python setup.py install

## Documentation

- **ADCP**
  - *adcp_binning()*: Binning ADCP data between CTD stations;
  - *mdr()*: Calibrate the normal geostrophic velocity section from observed normal velocity data;
- **AO**
  - *vectoa()*: Objective analysis for vectorial fields;
  - *scaloa()*: Objective analysis for scalar fields;
  - *Functions to calculate the correlation length*;
- **CTD**
  - *isopic_depth()*: Find isopycnal depths for density matrix;
  - *Processing functions*: Binning, despike, loopedit, hann filter, etc...;
- **DYN**
  - *dyn_amp()*: Makes the projection of every dynamical mode to velocity to obtain its amplitude;
  - *zeta()*: Calculates the vorticity field by velocity field;
  - *psi2uv()*: Calculates the velocity field by stream function scalar field;
  - *eqmodes()*: Calculates the equatorial dynamical pressure and vertical velocity modes; 
- **EOF**
  - *eoft()*: Calculates the Empirical Orthogonal Functions;
  - *my_eof_interp()*: Fillgaps on matrix based on EOFs (translated by Cesar Rocha Matlab version);
- **SEAPLOT**
  - *make_map()*: Create a map with inset axis for location (Based on Filipe Fernandes version);
- **utils**
  - *save_pickle()*: Save python object as pickle;
  - *load_pickle()*: Load python object by pickle file;
  - *interp2_yaxis()*: Interpolates a given data by second axis (Y-axis);
  - *rsme()*: Calculates the mean squared error;
  - *nans()*: Create an array full od NaNs with some given shape;
  - *argdistnear()*: Finds the index to nearest points in (xi,yi) from (x,y);
  - *download_bathy()*: Downloads ETOPO1 data and make subset;
   - *bathy_lims()*: Obtain bathymetry lines for boudary contition;
  - *deflagg()*: Remove badflags data from pandas.DataFrame;
  - *near()*: Finds the nearest value from list;
  - *argnear()*: Finds the index for nearest value from list;
  - *select_rad()*: Finds the index to select the section from 2 *matplotlib.pyplot.ginput()* points;
  - *extrap_all()*: Extrapolation function based on data's X gradient;

## Tests

Run the code below and it cannot return error:

python -c "from OceanLab import *"

## Contributors

Everyone can contribute to this code. Some of functions was based on Filipe Fernandes or Cesar Rocha functions and some of them was created with help of Hélio Almeida and Wandrey Watanabe at Laboratório de Dinâmica Oceânica of University of São Paulo (USP).

## License

MIT license.

##WARNING 
Be careful with the use of this. It can have errors.


Sao Paulo - SP, february 2016


[1]:http://iuryt.github.io

