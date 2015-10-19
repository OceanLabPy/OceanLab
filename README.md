# ocean
Package of Python scripts for Oceanography

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
		  
And needs some bathymetry and climatology data in particular format (mainly .mat from Matlab).
You can find it on the Example directory.

INSTALLATION:
I'm updating the package yet.

I will change it to get better but you have to clone the ocean directory and add to PYTHONPATH
the path to the directory that have it:

For example:
If ocean directory is in home folder (Linux). You can do it editing bashrc:

	oce=/home
	PYTHONPATH=$oce:$PYTHONPATH
	export PYTHONPATH
	
Sorry for the lack of comments on some functions but I'm working on it and this is unstable yet.
Until now, the less commented script is the seaplot.py

USAGE:

Some of directories has functions to import as: ADCP,AO,CTD,MOORING and SEAPLOT.

You can import as:

from ocean.AO import *

ans use like:

AO.psi2uv()

WARNING: Be careful with the use of this. It can have errors.

Sao Paulo - SP, march 2015


