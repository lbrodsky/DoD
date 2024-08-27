#!/usr/bin/env python3

# DEM (digital elevation model) of Difference (DoD)
# Python script for calculating glacier DoD indicators

import os
import glob
import argparse
import numpy as np
from osgeo import ogr
from osgeo import osr
from osgeo import gdal
gdal.UseExceptions()

__version__ = "0.1.0"

def volume_change():
    """...
    """

    return 0

def main(args):
    """...
    """

    pass



if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-y1', '--year1_dem', type=str, help='DEM for year 1.')
    parser.add_argument('-y2', '--year2_dem', type=str, help='DEM for year 1.')
    parser.add_argument('-d', '--dst_dir', metavar='dst_dir', type=str,
                        help='Target directory where the results are stored.')
    args = parser.parse_args()

    main(args)