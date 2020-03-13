"""
data.csatho : Csatho 2014 data import module.

This module provides functions to import surface elevation change (rates) data
from Csatho 2014 into a CISM dataset.

Functions list:
    *

Notes
-----
This data was downloaded from the citation below.

The data uses the EPSG:32624 (WGS84 / UTM zone 24N) projection:
    * Tranverse_Mercator
    * WGS84 ellipsoid
    * Latitude of projection origin = 0 degrees
    * Central meridian = -39 degrees
    * false eastings = 500000
    * false northings = 0
    * scale factor: 0.9996
    * Proj4: +proj=utm +zone=24 +datum=WGS84 +units=m +no_defs

Because the original data files are written in terms of ICESat Mission phases,
the data have been processed into a more obvious form. The processed data has
been distributed into analysis year NetCDF files, where the extension on the
.nc files indicate the years over which the rates of change were calculated
(e.g. "0304" were calculated between 2003 and 2004). Because the coordinates
appear to be transposed from standard convention in the original data, the
processed data's coordinates have been transposed (common features were
spot checked) and the thickness change variable has been renamed to dhdt.

References
----------
Beata M. Csatho, Anton F. Schenk, Cornelis J. van der Veen, Gregory Babonis,
Kyle Duncan, Soroush Rezvanbehbahani, Michiel R. van den Broeke, Sebastian B.
Simonsen, Sudhagar Nagarajan, Jan H. van Angelen: Greenland ice dynamics from
laser altimetry, Proceedings of the National Academy of Sciences Dec 2014,
111 (52) 18478-18483; DOI: 10.1073/pnas.1411680112

Data: https://arcticdata.io/catalog/view/doi:10.5065/D6HM56KS
"""

import numpy as np
import scipy
import scipy.interpolate
import pyproj
from util import speak
from util import projections
from util.ncfunc import copy_atts_add_fill
from build_antarctica import nn_interp

# ICESat campaign mid-dates in form of: {name: MMDDYY}
ICESAT_CAMPAIGN_MIDDATES = {
    "L1A": "030503",
    "L1B": "032503",
    "L2A": "100703",
    "L2B": "030504",
    "L2C": "060504",
    "L3A": "101904",
    "L3B": "030705",
    "L3C": "060605",
    "L3D": "110705",
    "L3E": "031106",
    "L3F": "060906",
    "L3G": "111006",
    "L3H": "032907",
    "L3I": "101907",
    "L3J": "030508",
    "L2D": "120608",
    "L2E": "031509",
    "L2F": "101509",
}

MISSING_VAL = np.float32(2.0e36)


def dhdt_epsg3413(args, nc_csatho, nc_base, base, proj_epsg3413):
    in_var = "L1A_to_L2F_Avg_dhdt"
    csatho = projections.DataGrid()
    csatho.lat = nc_csatho.variables["lat"]
    csatho.lon = nc_csatho.variables["lon"]
    csatho.dhdt = nc_csatho.variables[in_var]

    csatho.lat_grid, csatho.lon_grid = scipy.meshgrid(
        csatho.lat[:], csatho.lon[:], indexing="ij"
    )
    csatho.x, csatho.y = proj_epsg3413(
        csatho.lon_grid.ravel(), csatho.lat_grid.ravel()
    )

    speak.verbose(args, "    Interpolating dhdt and writing to base.")
    vy = csatho.y[~csatho.dhdt[:].mask.ravel()]
    vx = csatho.x[~csatho.dhdt[:].mask.ravel()]
    vz = csatho.dhdt[:].ravel()[~csatho.dhdt[:].mask.ravel()]

    dhdt = np.ma.masked_invalid(
        scipy.interpolate.griddata(
            np.column_stack((vy, vx)),
            vz,
            (base.y_grid, base.x_grid),
            method="linear",
        )
    )

    dhdt.mask = dhdt.mask | np.isclose(base.thk[:], 0.0)

    base.dhdt = nc_base.createVariable("dhdt", "f4", ("y", "x",))
    base.dhdt[:] = dhdt[:]
    # Copy the attributes, somehow it's missing the missing value setter
    copy_atts_add_fill(nc_csatho.variables[in_var], base.dhdt, MISSING_VAL)

    base.dhdt.long_name = "average {}--{} surface elevation change rate".format(
        ICESAT_CAMPAIGN_MIDDATES["L1A"][-2:],
        ICESAT_CAMPAIGN_MIDDATES["L2F"][-2:],
    )
    base.dhdt.standard_name = "tendency_of_land_ice_thickness"
    base.dhdt.units = "m year-1"
    base.dhdt.grid_mapping = "epsg_3413"
    base.dhdt.coordinates = "lon lat"
    base.dhdt.source = (
        "https://arcticdata.io/catalog/view/doi:10.5065/D6HM56KS"
    )
    base.dhdt.reference = (
        "Beata M. Csatho, Anton F. Schenk, Cornelis J. van der "
        "Veen, Gregory Babonis, Kyle Duncan, Soroush Rezvanbehbahani, "
        "Michiel R. van den Broeke, Sebastian B. Simonsen, Sudhagar Nagarajan, "
        "Jan H. van Angelen: Greenland ice dynamics from laser altimetry, "
        "Proceedings of the National Academy of Sciences Dec 2014, 111 (52) "
        "18478-18483; DOI: 10.1073/pnas.1411680112"
    )


def dhdt_bamber(
    args, nc_csatho, nc_base, base, proj_epsg3413, proj_eigen_gl04c
):
    in_var = "L1A_to_L2F_Avg_dhdt"
    csatho = projections.DataGrid()
    csatho.lat = nc_csatho.variables["lat"]
    csatho.lon = nc_csatho.variables["lon"]
    csatho.dhdt = nc_csatho.variables[in_var]

    csatho.lat_grid, csatho.lon_grid = np.meshgrid(
        csatho.lat[:], csatho.lon[:], indexing="ij"
    )

    csatho.x, csatho.y = proj_eigen_gl04c(
        csatho.lon_grid.ravel(), csatho.lat_grid.ravel()
    )

    speak.verbose(args, "    Interpolating dhdt and writing to base.")
    vy = csatho.y[~csatho.dhdt[:].mask.ravel()]
    vx = csatho.x[~csatho.dhdt[:].mask.ravel()]
    vz = csatho.dhdt[:].ravel()[~csatho.dhdt[:].mask.ravel()]

    dhdt = np.ma.masked_invalid(
        scipy.interpolate.griddata(
            np.column_stack((vy, vx)),
            vz,
            (base.y_grid, base.x_grid),
            method="linear",
        )
    )

    dhdt.mask |= np.isclose(base.thk[:], 0.0)

    dhdt = np.ma.masked_greater(dhdt, 1e10)
    dhdt = dhdt.filled(MISSING_VAL)
    base.dhdt = nc_base.createVariable("dhdt", "f4", ("y", "x",))
    base.dhdt[:] = dhdt[:]

    # Copy the attributes, somehow it's missing the missing value setter
    copy_atts_add_fill(nc_csatho.variables[in_var], base.dhdt, MISSING_VAL)

    # Now replace all the metadata except for the fill/missing value
    base.dhdt.long_name = "average {}--{} surface elevation change rate".format(
        ICESAT_CAMPAIGN_MIDDATES["L1A"][-2:],
        ICESAT_CAMPAIGN_MIDDATES["L2F"][-2:],
    )
    base.dhdt.standard_name = "tendency_of_land_ice_thickness"
    base.dhdt.units = "m year-1"
    base.dhdt.grid_mapping = "mcb"
    base.dhdt.coordinates = "lon lat"
    base.dhdt.source = (
        "https://arcticdata.io/catalog/view/doi:10.5065/D6HM56KS"
    )
    base.dhdt.reference = (
        "Beata M. Csatho, Anton F. Schenk, Cornelis J. van der "
        "Veen, Gregory Babonis, Kyle Duncan, Soroush Rezvanbehbahani, "
        "Michiel R. van den Broeke, Sebastian B. Simonsen, Sudhagar Nagarajan, "
        "Jan H. van Angelen: Greenland ice dynamics from laser altimetry, "
        "Proceedings of the National Academy of Sciences Dec 2014, 111 (52) "
        "18478-18483; DOI: 10.1073/pnas.1411680112"
    )


def dhdt_all(args, nc_csatho, nc_base, base, proj_out=None):
    """
    Interpolate dhdt data.

    Parameters
    ----------
    args :
    :param args:
    :param nc_csatho:
    :param nc_base:
    :param base:
    :param proj_out:
    :return:
    """
    prj_epsg, prj_mcb = projections.greenland()
    in_var = "dhdt"
    x_in, y_in = np.meshgrid(
        nc_csatho.variables["X"][:], nc_csatho.variables["Y"][:]
    )
    if proj_out is not None:
        # Transform first, data in is on Bamber (MCB) Grid.
        x_in, y_in = pyproj.transform(prj_mcb, proj_out, x=x_in, y=y_in)

    # dhdt = scipy.interpolate.griddata(
    #     np.column_stack((x_in.ravel(), y_in.ravel())),
    #     nc_csatho.variables[in_var][:].ravel(),
    #     (base.x_grid.ravel(), base.y_grid.ravel()),
    #     method="nearest",
    # )
    dhdt = nn_interp(
        x_in,
        y_in,
        base.x_grid,
        base.y_grid,
        nc_csatho.variables[in_var][:],
        nbrs=2,
    )
    dhdt = np.ma.masked_less(dhdt, -1e28)
    dhdt = dhdt.reshape(base.ny, base.nx)
    dhdt.mask |= np.isclose(base.thk[:], 0.0)

    dhdt = dhdt.filled(MISSING_VAL)

    base.dhdt = nc_base.createVariable("dhdt", "f4", ("y", "x",))
    base.dhdt[:] = dhdt[:]

    # Copy the attributes, somehow it's missing the missing value setter
    copy_atts_add_fill(nc_csatho.variables[in_var], base.dhdt, MISSING_VAL)

    # Now replace all the metadata except for the fill/missing value
    base.dhdt.long_name = "average {}--{} surface elevation change rate".format(
        ICESAT_CAMPAIGN_MIDDATES["L1A"][-2:],
        ICESAT_CAMPAIGN_MIDDATES["L2F"][-2:],
    )
    base.dhdt.standard_name = "tendency_of_land_ice_thickness"
    base.dhdt.units = "m year-1"
    base.dhdt.grid_mapping = "mcb"
    base.dhdt.coordinates = "lon lat"
    base.dhdt.source = (
        "https://arcticdata.io/catalog/view/doi:10.5065/D6HM56KS"
    )
    base.dhdt.reference = (
        "Beata M. Csatho, Anton F. Schenk, Cornelis J. van der "
        "Veen, Gregory Babonis, Kyle Duncan, Soroush Rezvanbehbahani, "
        "Michiel R. van den Broeke, Sebastian B. Simonsen, Sudhagar Nagarajan, "
        "Jan H. van Angelen: Greenland ice dynamics from laser altimetry, "
        "Proceedings of the National Academy of Sciences Dec 2014, 111 (52) "
        "18478-18483; DOI: 10.1073/pnas.1411680112"
    )
