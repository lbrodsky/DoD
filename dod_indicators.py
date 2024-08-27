#!/usr/bin/env python3

# Python script to calculate glacier DoD indicators

import os
import argparse
import numpy as np
from osgeo import gdal
gdal.UseExceptions()

import matplotlib.pyplot as plt

__version__ = "0.1.0"


def read_image(fn):
    """Read DoD to NumPy array.
    """
    ds = gdal.Open(fn)
    arr = ds.ReadAsArray()

    return arr

def area(data):
    """Calculate DoD area.
       (asuming data are in 1m pixel size)
       data = DoD
    """
    mask = np.isfinite(data)
    # plt.imshow(mask)
    count = np.sum(mask)
    print(f'{count} pixels (m^2)')

    return count


def volume_change(mean_change, area):
    """Total Volume Change [m3 x 106] = (Mean elevation change x Total area) / 106
    Let';s discuss with Susanne
    """
    return mean_change * area

def yearly_volume_change(total_volume_change, years):
    """Yearly Volume Change  [m3 x 106  yr-1] = Total Volume Change / Total years
    """
    return total_volume_change / years

def mass_change(total_volume_change, ICE_DENSITY=850):
    """Mass Change  [Kg] = (Volume change x 850)
    """
    mass_change = (total_volume_change * ICE_DENSITY)

    return mass_change

def mass_change_uncertainty(total_volume_change, ICE_DENSITY_UNCERTAINTY=60):
    """Mass Change uncertainty [Kg] = (Volume change x 850)
    """
    mass_change_uncertainty = (total_volume_change * ICE_DENSITY_UNCERTAINTY)

    return mass_change_uncertainty


def mass_change_unit(total_volume_change, ice_area, ICE_DENSITY=850):
    """Mass/unit area [kg/m2] = (Volume change x 850) / Total Area
    """
    mass_unit_area =  (total_volume_change * ICE_DENSITY) / ice_area

    return mass_unit_area

def mass_balance(mass_unit_area):
    """Mass balance [m w.e.] = Mass per unit area / 1000
    """
    mass_balance_ = mass_unit_area / 1000

    return mass_balance_

def annual_mass_change(mass_balance, years):
    """Annual mass change [m w.e. yr-1] = Mass balance / Total years
    """
    annual_mass_change_ = mass_balance / years

    return annual_mass_change_


def main(args):
    """...
    """

    print('Calculating DoD indicators!')
    print(f'Year 1 DEM: {os.path.basename(args.dod)}')
    print('---')

    DoD = read_image(args.dod)
    plt.hist(DoD)

    indicators = {}

    # 1. Mean elevation change
    mean_change = np.nanmean(DoD)
    indicators['mean_change'] = mean_change
    print(f'Mean change: {round(mean_change, 2)} m')

    median_change = np.nanmedian(DoD)
    indicators['median_change'] = median_change
    print(f'Median change: {round(median_change, 2)} m')

    # 2. Total Volume Change [m3 x 10^6]
    dod_area = area(DoD)
    total_volume_change = volume_change(mean_change, dod_area) / 10**6
    indicators['total_volume_change'] = total_volume_change
    print(f'Total Volume Change: {round(total_volume_change, 2)} [m3 x 10^6]')

    # 3. Yearly Volume Change  [m3 x 106  yr-1] = Total Volume Change / Total years
    volume_change_per_year = yearly_volume_change(total_volume_change, args.year2 - args.year1)
    indicators['yearly_volume_change'] = volume_change_per_year
    print(f'Yearly Volume Change: {round(volume_change_per_year, 2)} [m3 x 10^6 yr-1]')

    # 4. Mass Change  [Mt] = (Volume change x 850) / 109
    mass_change_mt = mass_change(total_volume_change, ICE_DENSITY=850) / 1000
    indicators['mass_change'] = mass_change_mt
    print(f'Mass Change: {round(mass_change_mt, 2)} [Mt]')

    # Uncertainty:
    mass_change_uncertainty_mt = mass_change_uncertainty(total_volume_change, ICE_DENSITY_UNCERTAINTY=60) / 1000
    indicators['mass_change_uncertainty'] = mass_change_uncertainty_mt
    print(f'Mass Change uncertainty: {round(mass_change_uncertainty_mt, 2)} [Mt]')

    # 5. Mass/unit area [kg/m2] = (Volume change x 850) / Total Area
    mass_unit_area = mass_change_unit(total_volume_change * 10**6, dod_area, ICE_DENSITY=850)
    indicators['mass_unit_area'] = mass_unit_area
    print(f'Mass / unit area: {round(mass_unit_area, 4)} [kg/m2]')

    # 6. Mass balance water equivalent [m w.e.] = Mass per unit area / 1000
    mass_balance_we = mass_balance(mass_unit_area)
    indicators['mass_balance'] = mass_balance_we
    print(f'Mass balance: {round(mass_balance_we, 4)} [m w.e.]')

    # 7. Annual mass change [m w.e. yr-1] = Mass balance / Total years
    annual_mass_change_we = annual_mass_change(mass_balance_we, (args.year2-args.year1))
    indicators['annual_mass_balance'] = annual_mass_change_we
    print(f'Annual Mass balance: {round(annual_mass_change_we, 4)} [m w.e. yr-1]')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dod', type=str, help='DEM of Differnce.')
    parser.add_argument('-i', '--ice_density', type=int, help='Ice density constant in kg/m^3.',
                        default=850)
    parser.add_argument('-y1', '--year1', type=int, help='Year 1 DEM.')
    parser.add_argument('-y2', '--year2', type=int, help='Year 2 DEM.')
    args = parser.parse_args()

    """
    class Namespace:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    args = Namespace(dod='/Users/lukas/Work/prfuk/clanky/AUC_Special_issue/Belvedere_3D_change/data/DSMs_1m/DoD/Belvedere_2009_1951_DoD.tif', 
        year1 = 1951,
        year2 = 2009 
    ) 
    """

    main(args)