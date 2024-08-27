#!/usr/bin/env python3

# Python script for DEM (Digital Elevation Model) of Difference (DoD) calculation
# Assumption: DEMs (DSMs) are cor-registered and use the same extent and pixel size
# Example usage:
"""
python3 dod_calculation.py \
    -y1 ./data/DSMs_1m/Belvedere_1951_DEM_1m_ext09.tif \
    -y2 ./data/DSMs_1m/Belvedere_2009_DEM_1m_ext09.tif \
    -d ./data/DSMs_1m/DoD
"""
# TODO: create procedure to spatially align the two DEMs (DSMs) before DoD


import os
import argparse
import numpy as np
from osgeo import gdal
from osgeo import gdal_array
import matplotlib.pylab as plt

gdal.UseExceptions()

__version__ = "0.1.0"


def read_image(fn):
    """Read image/DEM to NumPy array.
    """
    ds = gdal.Open(fn)
    dem_arr = ds.ReadAsArray()

    return dem_arr

def dod(dem1, dem2):
    """Simple DEM difference
    """
    return dem1 - dem2

def mask_dem(dem, mask_fn):
    """Masking DEM
    """
    mask = read_image(mask_fn)

    return dem * mask

def none_treat(dem, ndvs=[-0, -99]):
    """None value tretment.
       No data value: 0
    """
    for val in ndvs:
        # print(val)
        np.place(dem, dem==val, [0])

    dem[dem == 0] = np.nan

    return dem

def median_filter(data, mask):
    """Median 2D filter for masked values.
    """
    print('Median filter started.')
    m, n = data.shape
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            if mask[i, j]:
                # So far only 3x3 kernel
                temp = [data[i - 1, j - 1],
                        data[i - 1, j],
                        data[i - 1, j + 1],
                        data[i, j - 1],
                        data[i, j],
                        data[i, j + 1],
                        data[i + 1, j - 1],
                        data[i + 1, j],
                        data[i + 1, j + 1]]
                temp = sorted(temp)
                data[i, j] = temp[4]
    print('Finished.')

    return data

def write_imge(DoD, args):
    """Write resulting DoD image based inputs.
    """
    tokens1 = os.path.basename(args.year1_dem).split('_')
    tokens2 = os.path.basename(args.year2_dem).split('_')
    file_name = os.path.join(args.dst_dir, '_'.join([tokens1[0], tokens2[1], tokens1[1], 'DoD']) + '.tif')

    # geo metadata
    ds = gdal.Open(args.year1_dem, gdal.GA_ReadOnly)
    projection = ds.GetProjection()
    geo_transform = ds.GetGeoTransform()
    nrow = ds.RasterYSize
    ncol = ds.RasterXSize
    nband = ds.RasterCount
    gdal_frmt = "GTiff"
    dtype = np.float64
    gdal_dtype = gdal_array.NumericTypeCodeToGDALTypeCode(dtype)
    ndv = -999
    driver = gdal.GetDriverByName(str(gdal_frmt))

    ds_out = driver.Create(file_name, ncol, nrow, 1, gdal_dtype)
    ds_out.GetRasterBand(1).WriteArray(DoD)
    ds_out.GetRasterBand(1).SetNoDataValue(ndv)
    ds_out.SetProjection(projection)
    ds_out.SetGeoTransform(geo_transform)
    ds_out = None

    print(f'File written: {os.path.basename(file_name)}')


def main(args):
    """Man DoD process.
    """
    print('Calculating DoD!')
    print(f'Year 1 DEM: {os.path.basename(args.year1_dem)}')
    print(f'Year 2 DEM: {os.path.basename(args.year2_dem)}')
    print('---')

    # Read DEMs
    y1_dem = read_image(args.year1_dem)
    y2_dem = read_image(args.year2_dem)
    # Masks
    y1_dem_msk = mask_dem(y1_dem, args.mask1)
    y2_dem_msk = mask_dem(y2_dem, args.mask1)

    # Align Nones
    y1_dem_msk_none = none_treat(y1_dem_msk, ndvs=[-0, -99])
    y2_dem_msk_none = none_treat(y2_dem_msk, ndvs=[-0, -99])

    # Calculate DoD
    DoD = dod(y2_dem_msk_none, y1_dem_msk_none)

    # DoD confidence intervals
    sample_mean = np.nanmean(DoD)
    sigma = np.nanstd(DoD)
    c_low = sample_mean - (3 * sigma)
    c_high = sample_mean + (3 * sigma)

    # Plot image difference
    # max_int = max(abs(c_low), c_high)
    # plt.imshow(DoD, vmin=-max_int/1.5, vmax=max_int/1.5, cmap='RdBu')
    # plt.colorbar()

    # Filter outliers
    msk = np.logical_and(np.logical_or(DoD < c_low, DoD > c_high), np.isfinite(DoD))
    DoD_fil = median_filter(DoD, msk)

    # Plot
    # plt.imshow(DoD_fil, vmin=-max_int/1.5, vmax=max_int/1.5, cmap='RdBu')
    # plt.colorbar()
    print(f'Median change: {np.nanmedian(DoD_fil)}')

    # Save DoD to file
    write_imge(DoD_fil, args)
    print('---')


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-y1', '--year1_dem', type=str, help='DEM for year 1.')
    parser.add_argument('-y2', '--year2_dem', type=str, help='DEM for year 2.')
    parser.add_argument('-m1', '--mask1', type=str, help='Mask of glacier for year 1.')
    parser.add_argument('-m2', '--mask2', type=str, help='Mask of glacier for year 2.')
    parser.add_argument('-d', '--dst_dir', metavar='dst_dir', type=str,
                        help='Target directory where the results are stored.')
    parser.add_argument('-t', '--year_token', type=int, default=1, help='Index of year token.')
    args = parser.parse_args()

    """
    class Namespace:
        def __init__(self, **kwargs):
            self.__dict__.update(kwargs)

    args = Namespace(year1_dem='/Users/lukas/Work/prfuk/clanky/AUC_Special_issue/Belvedere_3D_change/data/DSMs_1m/Belvedere_1951_DEM_1m_ext09.tif',
                     year2_dem='/Users/lukas/Work/prfuk/clanky/AUC_Special_issue/Belvedere_3D_change/data/DSMs_1m/Belvedere_2009_DEM_1m_ext09.tif',
                     mask1='/Users/lukas/Work/prfuk/clanky/AUC_Special_issue/Belvedere_3D_change/data/DSMs_1m/Belvedere_1951_mask_1m_ext09.tif',
                     mask2='/Users/lukas/Work/prfuk/clanky/AUC_Special_issue/Belvedere_3D_change/data/DSMs_1m/Belvedere_2015_mask_1m_ext09.tif',
                     dst_dir='/Users/lukas/Work/prfuk/clanky/AUC_Special_issue/Belvedere_3D_change/data/DSMs_1m/DoD')
    """

    main(args)
