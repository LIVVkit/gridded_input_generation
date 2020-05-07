"""
data.racmo2p3 : RACMO 2.3 data import module.

This module provides functions to import data from a 1km downscaled RACMO 2.3 
dataset into a CISM dataset. 

Functions list:
    * acab_epsg3413(args, nc_racmo, nc_base, base)

Notes
-----
This data was provided by Jan Lenaerts from the citation below.  

The data uses the ESPG:3413 projection:
    * Polar stereographic
    * WGS84 ellipsoid
    * Standard parallel = 70 degrees
    * Latitude of projection origin = 90 degrees
    * Central meridian = -45 degrees
    * false eastings = 0
    * false northings = 0


References
----------
Noel, B., van de Berg, W. J., Machguth, H., Lhermitte, S., Howat, I., Fettweis,
X., and van den Broeke, M. R.: a daily, 1 km resolution data set of downscaled
Greenland ice sheet surface mass balance (1958--2015), The Cryosphere, 10,
2361-2377, doi:10.5194/tc-10-2361-2016, 2016. 
"""

import sys
import scipy
import numpy as np
import pyproj

from util.ncfunc import copy_atts_bad_fill
from util import speak
from util import projections
from pykdtree.kdtree import KDTree

MISSING_VAL = np.float32(2.0e36)


def acab_epsg3413(args, nc_racmo, nc_base, base):
    racmo = projections.DataGrid()
    racmo.y = nc_racmo.variables["y"]
    racmo.x = nc_racmo.variables["x"]
    racmo.ny = racmo.y[:].shape[0]
    racmo.nx = racmo.x[:].shape[0]
    racmo.make_grid()

    base_vars = {"smb": "SMB_mean", "smb_std": "SMB_stdv"}
    for bvar, rvar in base_vars.items():
        speak.verbose(
            args, "   Interpolating " + bvar + " and writing to base."
        )
        sys.stdout.write("   [%-60s] %d%%" % ("=" * 0, 0.0))
        sys.stdout.flush()
        racmo_data = np.ma.masked_values(
            nc_racmo.variables[rvar][-1, :, :], 9.96921e36
        )
        racmo_data *= 12  # Convert from (mm w.e. / mon) to (mm w.e. / year)

        data_min = racmo_data.min()
        data_max = racmo_data.max()

        racmo_to_base = scipy.interpolate.RectBivariateSpline(
            racmo.y[:], racmo.x[:], racmo_data[:, :], kx=1, ky=1, s=0
        )  # regular 2d linear interp. but faster
        base_data = np.zeros(base.dims)
        for ii in range(0, base.nx):
            ctr = (ii * 60) // base.nx
            sys.stdout.write(
                "\r   [%-60s] %d%%" % ("=" * ctr, ctr / 60.0 * 100.0)
            )
            sys.stdout.flush()
            base_data[:, ii] = racmo_to_base.ev(
                base.y_grid[:, ii], base.x_grid[:, ii]
            )
        sys.stdout.write("\r   [%-60s] %d%%\n" % ("=" * 60, 100.0))
        sys.stdout.flush()

        base_data[base_data < data_min] = MISSING_VAL
        base_data[base_data > data_max] = MISSING_VAL

        base.var = nc_base.createVariable(bvar, "f4", ("y", "x",))
        base.var[:] = base_data[:]
        copy_atts_bad_fill(nc_racmo.variables[rvar], base.var, MISSING_VAL)
        base.var.long_name = "Water Equivalent Surface Mass Balance"
        base.var.standard_name = "land_ice_lwe_surface_specific_mass_balance"
        if "std" in bvar:
            base.var.long_name += " standard deviation"
            base.var.standard_name += " standard_deviation"
        base.var.units = "mm year-1"
        base.var.grid_mapping = "epsg_3413"
        base.var.coordinates = "lon lat"
        base.var.source = "Brice Nöel"
        base.var.reference = (
            "Noël, Brice and van de Berg, Willem Jan and Lhermitte, Stef and "
            "van den Broeke, Michiel R. 2019: Rapid ablation zone expansion "
            "amplifies north Greenland mass loss. Science Advances."
            "doi:10.1126/sciadv.aaw0123"
        )


def acab_bamber(args, nc_racmo2p3, nc_base, base):
    """Get acab from the RACMO 2.0 data.

    This function pulls in the `smb` variable from the RACMO 2.0 dataset
    and writes it to the base dataset as `acab`. NetCDF attributes are
    preserved.
    Parameters
    ----------
    args :
        Namespace() object holding parsed command line arguments.
    nc_racmo2p3 :
        An opened netCDF Dataset containing the RACMO 2.0 data.
    nc_base :
        The created netCDF Dataset that will contain the base data.
    base :
        A DataGrid() class instance that holds the base data grid information.

    """
    x_in = nc_racmo2p3.variables["x"][:]
    y_in = nc_racmo2p3.variables["y"][:]
    x2d, y2d = np.meshgrid(x_in, y_in)

    proj_greenland = projections.greenland()
    base_vars = {"smb": "SMB_mean", "smb_std": "SMB_stdv"}

    x_t, y_t = pyproj.transform(
        proj_greenland[0], proj_greenland[1], x=x2d.flatten(), y=y2d.flatten(),
    )

    # transform_tree = scipy.spatial.cKDTree(np.vstack((x_t, y_t)).T)
    transform_tree = KDTree(np.vstack((x_t, y_t)).T)

    qd, qi = transform_tree.query(
        np.vstack((base.x_grid.flatten(), base.y_grid.flatten())).T, k=1,
    )

    for bvar, rvar in base_vars.items():
        speak.verbose(
            args, f"   Interpolating {rvar} and writing {bvar} to base."
        )
        z_in = nc_racmo2p3.variables[rvar][0, ...]
        z_in *= 12  # Convert from (mm w.e. / mon) to (mm w.e. / year)

        x_t, y_t, z_t = pyproj.transform(
            proj_greenland[0],
            proj_greenland[1],
            x=x2d.flatten(),
            y=y2d.flatten(),
            z=z_in.flatten(),
        )

        z_t = z_t.reshape(z_in.shape)
        z_t = z_t.filled(MISSING_VAL)

        z_interp = z_t.flatten()[qi].reshape(base.y_grid.shape)
        base.var = nc_base.createVariable(bvar, "f4", ("y", "x",))
        base.var[:] = z_interp[:]

        # If there are missing values from the original dataset, update the
        # fill value to something easy
        try:
            _mvin = nc_racmo2p3.variables[rvar].getncattr("_FillValue")
        except AttributeError:
            _mvin = MISSING_VAL
        z_interp[np.where(z_interp == _mvin)] = MISSING_VAL

        copy_atts_bad_fill(nc_racmo2p3.variables[rvar], base.var, MISSING_VAL)

        base.var.long_name = "Water Equivalent Surface Mass Balance"
        base.var.standard_name = "land_ice_lwe_surface_specific_mass_balance"
        if "std" in bvar:
            base.var.long_name += " standard deviation"
            base.var.standard_name += " standard_deviation"

        base.var.units = "mm year-1"
        base.var.grid_mapping = "mcb"
        base.var.coordinates = "lon lat"
        base.var.source = "Brice Nöel"
        base.var.reference = (
            "Noël, Brice and van de Berg, Willem Jan and Lhermitte, Stef and "
            "van den Broeke, Michiel R. 2019: Rapid ablation zone expansion "
            "amplifies north Greenland mass loss. Science Advances."
            "doi:10.1126/sciadv.aaw0123"
        )
        base.var.comments = "Nearest neighbour remapping"
