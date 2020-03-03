# -*- coding: utf-8 -*-
"""
MEASURES velocity data data.
Data is on EPSG:3414
https://nsidc.org/data/NSIDC-0725/versions/1

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
import pyproj
from pykdtree.kdtree import KDTree

from util import speak
from util import projections

DATA_PATH = "data"
MISSING_VAL_INPUT = -2e9
MISSING_VAL_OUT = np.float32(2e36)

META_ALL = {
    "source": "https://nsidc.org/data/NSIDC-0725/versions/1",
    "reference": (
        "Joughin, I. 2017, updated 2019. MEaSUREs Greenland Annual Ice "
        "Sheet Velocity Mosaics from SAR and Landsat, Version 1. "
        "Boulder, Colorado USA. NASA National Snow and Ice Data Center "
        "Distributed Active Archive Center. "
        "doi: https://doi.org/10.5067/OBXCG75U7540."
    ),
    "units": "m year-1",
    "missing_value": MISSING_VAL_OUT,
}

META_VAR = {
    "vx": {
        "long_name": "surface x velocity",
        "standard_name": "land_ice_surface_x_velocity",
    },
    "vy": {
        "long_name": "surface y velocity",
        "standard_name": "land_ice_surface_y_velocity",
    },
    "ex": {
        "long_name": "surface x velocity error",
        "standard_name": "land_ice_surface_x_velocity standard_error",
    },
    "ey": {
        "long_name": "surface y velocity error",
        "standard_name": "land_ice_surface_y_velocity standard_error",
    },
}


def velocity(args, nc_base, base, new_proj=None):
    """
    Interpolate (and transform if needed) Measures velocity data to 1km grid.

    Parameters
    ----------
    args : argparser.arguments
        Command line arguments
    nc_base : netCDF4.Dataset
        Opened dataset
    base : projections.DataGrid
        Output data grid specification
    new_proj : proj4.projection
        Projection of `base` grid, if different from EPSG:3413

    """
    vel_vars = ["vx", "vy", "ex", "ey"]
    _inxvar = "x"
    _inyvar = "y"
    prj_epsg_3413, _ = projections.greenland()

    if new_proj is None:
        META_ALL["grid_mapping"] = "epsg_3413"
    else:
        META_ALL["grid_mapping"] = "mcb"

    for vel_var in vel_vars:
        in_file = Path(
            DATA_PATH,
            "MEASURES_annual_vel_mosaics_v01",
            f"greenland_vel_mosaic500_2018_{vel_var}_v01.tif",
        )
        speak.verbose(args, f"   opened {vel_var}: {in_file}")

        data_in = xr.open_rasterio(in_file)
        in_x, in_y = np.meshgrid(data_in[_inxvar], data_in[_inyvar])
        data_in = data_in.where(data_in != MISSING_VAL_INPUT).values
        if new_proj is not None:
            speak.verbose(args, "     transform data to Bamber grid")
            in_x, in_y, data_in = pyproj.transform(
                prj_epsg_3413,
                new_proj,
                x=in_x.ravel(),
                y=in_y.ravel(),
                z=data_in.flatten(),
            )

        tree = KDTree(np.vstack([in_x.ravel(), in_y.ravel()]).T)

        _, idx = tree.query(
            np.vstack([base.x_grid.ravel(), base.y_grid.ravel()]).T
        )
        data_out = data_in.ravel()[idx].reshape(base.ny, base.nx)
        ncvar = nc_base.createVariable(vel_var, "f4", ("y", "x",))
        ncvar[:] = data_out
        ncvar.setncatts(META_ALL)
        ncvar.setncatts(META_VAR[vel_var])
