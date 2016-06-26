# -*- coding: utf-8 -*-


__version__ = '2.0'
__doc__='''This lib was made by Iury Sousa (simoesiury@gmail.com) from Brazil to help him and other people that
		  needs to organize, process and plot CTD and VADCP data from cruise. Manly directioned
		  to mesoscale oceanographers. The scripts need to have installed some packages:
		  
		  Matplotlib
		  Basemap
		  Numpy
		  Scipy
		  Pandas
		  Gibbs SeaWater
		  SeaWater
		  
		  And needs some bathymetry and climatology data in particular format (mainly .mat from Matlab)
		  Sao Paulo - SP, march 2015'''

from seaext import *
from seacalc import *
from seaplot import *
from misc import *