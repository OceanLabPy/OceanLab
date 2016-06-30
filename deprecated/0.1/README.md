# OceanLab
Package of Python scripts for Oceanography - Version 0.1

This lib was made by Iury Sousa (simoesiury@gmail.com) from Brazil to help him and other people that
needs to organize, process and plot CTD and VADCP data from cruise. Manly directioned
to mesoscale oceanographers. The scripts need to have installed some packages:
		  
		  Matplotlib
		  Basemap
		  Numpy
		  Scipy
		  Pandas
		  Gibbs SeaWater
		  SeaWater
		  
		  
INSTALLATION:
I'm updating the package yet.

Clone this directory and run:

python setup.py install
	
Sorry for the lack of comments on some functions but I'm working on it and this is unstable yet.
Until now, the less commented script is the SEAPLOT.py

USAGE:

Some of directories has functions to import as: ADCP,AO,CTD,EOF_fillgap,MOORING and SEAPLOT.

You can import as:

from OceanLab import AO

and use as:

AO.psi2uv()


WARNING: Be careful with the use of this. It can have errors.
Sao Paulo - SP, february 2016


