# -*- coding: utf-8 -*-
"""
MEASURES velocity data data.
Data is on EPSG:3414
https://nsidc.org/data/NSIDC-0725/versions/1?qt-data_set_tabs=1#qt-data_set_tabs

This data set, part of the NASA Making Earth System Data Records for Use in
Research Environments (MEaSUREs) program, contains annual ice velocity mosaics
for the Greenland Ice Sheet derived from Synthetic Aperture Radar (SAR) data
obtained by the German Space Agency's TerraSAR-X/TanDEM-X (TSX/TDX) and the
European Space Agency's Copernicus Sentinel-1A and -1B satellites, and from
the US Geological Survey's Landsat 8 optical imagery for years 2015 to 2018.

Data are provided at both 200 m and 500 m postings in GeoTIFF (.tif) format.
Six data files are available for each data year and resolution: a velocity
magnitude map (vv); separate x- and y-component velocities (vx, vy);
separate x- and y-component error estimates (ex, ey); and a browse file
(with color log scale velocity saturating at 3000 m/year). In addition,
ancillary files are provided for each data year to indicate the source
image pairs that were processed to produce the mosaics. These are provided
in two separate shapefiles (.shp) for the US Geological Survey (USGS)-provided
Landsat 8 (L8) and the German Aerospace Center (DLR) and European Space
Agency (ESA)-provided Synthetic Aperture Radar (SAR) data.

"""
from pathlib import Path
import numpy as np
import scipy
import scipy.interpolate
import xarray as xr

from util import speak
from util import projections


def velocity_epsg3413(args, base, proj_epsg3413):
    """

    Parameters
    ----------
    args : argparser.arguments
        Command line arguments
    base : projections.DataGrid
        Output data grid specification
    proj_epsg3413 : proj4.projection
        EPSG:3413 projection

    """
    vars = ["vx", "vy", "ex", "ey"]
    for vel_var in vars:
        in_file = Path(
            "MEASURES_annual_vel_mosaics_v01",
            f"greenland_vel_mosaic500_2018_{vel_var}_v01.tif",
        )
        data_in = xr.open_rasterio(in_file)

        # Data needs to be masked, then sampled at the new grid
        # no transform required since it's already on the EPSG:3413 grid
