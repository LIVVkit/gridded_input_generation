#!/usr/bin/env python2
# encoding: utf-8

"""
Build the antarctica dataset quick and dirty like based on an old version
"""
import os
import sys

import pyproj
import scipy.interpolate
import numpy as np

from argparse import Namespace
from datetime import datetime

from nco import Nco

from util import ncfunc
from util import projections
from util import finalize
import xarray as xr


FILL = True


def bedmachine_interp(x, y, data, new):
    """Interpolate gridded data on rectangular grid with Bivariate spline."""
    if y[1] - y[0] < 0:
        _yslice = slice(None, None, -1)
    else:
        _yslice = slice(None, None, None)

    if x[1] - x[0] < 0:
        _xslice = slice(None, None, -1)
    else:
        _xslice = slice(None, None, None)

    to_new = scipy.interpolate.RectBivariateSpline(
        y[_yslice], x[_xslice], data[_yslice, _xslice], kx=1, ky=1, s=0
    )
    new_data = np.zeros(new.dims)

    for ii in range(0, new.nx):
        new_data[:, ii] = to_new.ev(new.y_grid[:, ii], new.x_grid[:, ii])

    # new_data = new_data[_yslice, _xslice]
    return new_data


def main(coarse_out=False):
    """Check data files, load them, interpolate, output CISM input data."""
    # Data files Antarctica
    f_1km_out = (
        "ncs/antarctica_1km_" + datetime.now().strftime("%Y_%m_%d") + ".nc"
    )
    f_1km_in = "ncs/antarctica_1km_2017_05_03.nc"
    subshelf = (
        "data/Rignot-subShelfMeltRates/Ant_MeltingRate.flipNY.newAxes.nc"
    )

    bed_opt = "BedMachine"
    if bed_opt == "BedMachine":
        bedmap = (
            "data/500m.MassConsBed.AIS.Morlighem.2019/"
            "BedMachineAntarctica_2019-11-05_v01.nc"
        )
        bedmap_uncert = bedmap

    else:
        bedmap = "data/BISICLES/Antarctica-1km.BISICLES.CISM-style.nc"
        bedmap_uncert = (
            "data/5km-res-Bedmap2-ThkError/"
            "Antarctica-Bedmap2-thicknessUncertainty.nc"
        )

    heat_flux_file = "data/Martos-AIS-heatFlux/Antarctic_GHF.xyz"
    heat_flux_unc_file = (
        "data/Martos-AIS-heatFlux/Antarctic_GHF_uncertainty.xyz"
    )
    cryosat_file = "data/Cryosat2/CS2_dzdt.nc"
    all_found = True

    input_files = [
        f_1km_in,
        subshelf,
        bedmap,
        bedmap_uncert,
        heat_flux_file,
        heat_flux_unc_file,
        cryosat_file,
        "util/egm08_25.gtx",
    ]

    for _file in input_files:
        if not os.path.exists(_file):
            print("FILE NOT FOUND: " + _file)
            all_found = False
        else:
            print("  FILE FOUND: " + _file)

    if not all_found:
        sys.exit(1)

    f_template = ""
    coarse_list = [2, 4, 5, 8]
    args = Namespace()
    args.quite = True
    args.verbose = False

    # import shutil
    # shutil.copy('antarctica_1km_2017_05_03.nc', f_1km)
    print("NCKS -> 1 km")
    nco = Nco()
    nco.ncks(
        input=f_1km_in,
        output=f_1km_out,
        options=["-x", "-v", "acab_alb,artm_alb,dzdt"],
    )
    print("\nmake new grid " + f_1km_out)
    nc_new = ncfunc.get_nc_file(f_1km_out, "r+")
    new = projections.DataGrid()
    new.y = nc_new.variables["y1"]
    new.x = nc_new.variables["x1"]
    new.ny = new.y.shape[0]
    new.nx = new.x.shape[0]
    new.make_grid()

    print("\nsubshelf melt rates")
    nc_subm = ncfunc.get_nc_file(subshelf, "r")
    subm = projections.DataGrid()
    subm.y = nc_subm.variables["y1"]
    subm.x = nc_subm.variables["x1"]
    subm.ny = subm.y[:].shape[0]
    subm.nx = subm.x[:].shape[0]
    subm.make_grid()
    print("  Interpolate subm")
    subm_actu_to_new = scipy.interpolate.RectBivariateSpline(
        subm.y[:],
        subm.x[:],
        nc_subm.variables["melt_actual"][:, :],
        kx=1,
        ky=1,
        s=0,
    )
    new_subm_actu = np.zeros(new.dims)
    for ii in range(0, new.nx):
        new_subm_actu[:, ii] = subm_actu_to_new.ev(
            new.y_grid[:, ii], new.x_grid[:, ii]
        )

    nc_new.createVariable("subm", "f4", ("time", "y1", "x1",))
    nc_new.variables["subm"][0, :, :] = new_subm_actu[:, :]
    nc_new.variables["subm"].long_name = "sub-shelf melt rate"
    nc_new.variables["subm"].units = "m year-1"
    nc_new.variables["subm"].grid_mapping = "epsg_3031"
    nc_new.variables["subm"].coordinates = "lon lat"
    nc_new.variables["subm"].reference = (
        "Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, 2013: Ice-Shelf "
        "Melting Around Antarctica. Science, 341, 266-270, "
        "doi:10.1126/science.1235798."
    )
    nc_new.variables["subm"].source = "J. Mouginot"
    nc_new.variables["subm"].comments = (
        "2D linear interpolation of provided dataset. Authors request we "
        "communicate with them prior to publishing and acknowledge them as the "
        "source of the data in presentations."
    )
    # FIXME nc_new.variables['subm'].missing_value = FIXME

    subm_stdy_to_new = scipy.interpolate.RectBivariateSpline(
        subm.y[:],
        subm.x[:],
        nc_subm.variables["melt_steadystate"][:, :],
        kx=1,
        ky=1,
        s=0,
    )
    new_subm_stdy = np.zeros(new.dims)
    for ii in range(0, new.nx):
        new_subm_stdy[:, ii] = subm_stdy_to_new.ev(
            new.y_grid[:, ii], new.x_grid[:, ii]
        )

    nc_new.createVariable("subm_ss", "f4", ("time", "y1", "x1",))
    nc_new.variables["subm_ss"][0, :, :] = new_subm_stdy[:, :]
    nc_new.variables["subm_ss"].long_name = "steady state sub-shelf melt rate"
    nc_new.variables["subm_ss"].units = "m year-1"
    nc_new.variables["subm_ss"].grid_mapping = "epsg_3031"
    nc_new.variables["subm_ss"].coordinates = "lon lat"
    nc_new.variables["subm_ss"].reference = (
        "Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, 2013: Ice-Shelf "
        "Melting Around Antarctica. Science, 341, 266-270, "
        "doi:10.1126/science.1235798."
    )
    nc_new.variables["subm_ss"].source = "J. Mouginot"
    nc_new.variables["subm_ss"].comments = (
        "2D linear interpolation of provided dataset. Authors request we "
        "communicate with them prior to publishing and acknowledge them as the "
        "source of the data in presentations."
    )
    # FIXME nc_new.variables['subm_ss'].missing_value = FIXME
    print("\nBedmap2")
    # FIXME: Do the cleanup myself from bedmap2 data.
    _bedx = "x"
    _bedy = "y"
    _thck = "thickness"

    # nc_bike = ncfunc.get_nc_file(bedmap, "r")
    # bike = projections.DataGrid()
    # bike.y = nc_bike.variables[_bedy]
    # bike.x = nc_bike.variables[_bedx]
    # bike.ny = bike.y[:].shape[0]
    # bike.nx = bike.x[:].shape[0]
    # bike.make_grid()

    # print("  Interpolate thk")
    # new_bike = bedmachine_interp(bike.x, bike.y, nc_bike.variables[_thck], new)
    nc_bike = xr.open_dataset(bedmap)
    new_bike = nc_bike[_thck].interp(
        **{_bedx: new.x[:], _bedy: new.y[:]}, method="nearest"
    )
    new_bike = new_bike.values

    if FILL:
        new_bike = np.ma.masked_where(new_bike < 11.0, new_bike)
        nc_new.variables["thk"][0, :, :] = new_bike[:, :].filled(0.0)
    else:
        nc_new.variables["thk"][0, :, :] = new_bike[:, :]

    if bed_opt == "BedMachine":
        meta = {
            "source": (nc_bike.Title + " " + nc_bike.Author),
            "reference": (
                nc_bike.Data_citation + ", " + "doi:10.1038/s41561-019-0510-8"
            ),
            "comments": (
                "Obtained from NSIDC: https://nsidc.org/data/nsidc-0756. "
                "Resampled from 500m grid using BivariateSpline to 1km."
            ),
        }

    else:
        meta = {
            "source": (
                "Antarctica-1km.BISICLES.CISM-style.nc; generated by "
                "Dan Martin or maybe Steph Cornford)."
            ),
            "reference": (
                "Fretwell, P., et al.: Bedmap2: improved ice bed, "
                "surface and thickness datasets for Antarctica, The Cryosphere,"
                " 7, 375-393, https://doi.org/10.5194/tc-7-375-2013, 2013."
            ),
            "comments": (
                "Origionally from the BEedmap2 dataset, but the water portion "
                "of Lake Vostok has been filled with ice and the peripheral "
                "islands have been removed so ice thickness only represents "
                "the main ice sheet."
            ),
        }

    nc_new.variables["thk"].long_name = "ice thickness"
    nc_new.variables["thk"].standard_name = "land_ice_thickness"
    nc_new.variables["thk"].units = "m"
    nc_new.variables["thk"].ancillary_variables = "thkerr"
    nc_new.variables["thk"].grid_mapping = "epsg_3031"
    nc_new.variables["thk"].coordinates = "lon lat"
    for meta_var in meta:
        nc_new.variables["thk"].setncattr(meta_var, meta[meta_var])

    # FIXME nc_new.variables['thk'].missing_value = FIXME
    print("\nBedmap2 Uncertainty")
    nc_berr = ncfunc.get_nc_file(bedmap_uncert, "r")
    _berrvar = "errbed"
    berr = projections.DataGrid()
    berr.x = nc_berr.variables[_bedx][:]
    berr.y = nc_berr.variables[_bedy][:]
    berr.ny = berr.y.shape[0]
    berr.nx = berr.x.shape[0]
    berr.make_grid()

    print("  Interpolate thkerr")
    # new_berr = bedmachine_interp(berr.y, berr.x, nc_berr[_berrvar], new)
    new_berr = nc_bike[_berrvar].interp(
        **{_bedx: new.x[:], _bedy: new.y[:]}, method="nearest"
    )
    new_berr = new_berr.values

    # Mask out below 0s
    new_berr = np.ma.masked_less(new_berr, 0)

    if bedmap_uncert != bedmap:
        berr.points = list(zip(berr.y_grid.flatten(), berr.x_grid.flatten()))
        berr_to_new = scipy.interpolate.NearestNDInterpolator(
            berr.points, nc_berr[_berrvar][:].flatten()
        )
        new_berr = np.zeros(new.dims)
        for ii in range(0, new.nx):
            new_berr[:, ii] = berr_to_new(
                list(zip(new.y_grid[:, ii], new.x_grid[:, ii]))
            )

    nc_new.createVariable("thkerr", "f4", ("time", "y1", "x1",))
    if FILL:
        new_berr = np.ma.masked_values(new_berr, -9999.0)
        nc_new.variables["thkerr"][0, :, :] = new_berr[:, :].filled(150.0)
    else:
        nc_new.variables["thkerr"][0, :, :] = new_berr[:, :]

    nc_new.variables["thkerr"].long_name = "ice thickness error"
    nc_new.variables[
        "thkerr"
    ].standard_name = "land_ice_thickness standard_error"
    nc_new.variables["thkerr"].units = "m"
    nc_new.variables["thkerr"].grid_mapping = "epsg_3031"
    nc_new.variables["thkerr"].coordinates = "lon lat"

    nc_new.variables["thkerr"].reference = meta["reference"]
    if bed_opt == "BedMachine":
        nc_new.variables[
            "thkerr"
        ].comments = "Resampled from 500m grid using BivariateSpline to 1km."
    else:
        nc_new.variables["thkerr"].comments = (
            "Resampled from 5km grid using 2D nearest neighbor interpolation, "
            "missing values filled with a value of 150.0. "
        )

    if not FILL:
        nc_new.variables["thkerr"].missing_value = -9999.0

    print("\nHeat Flux")
    heat_flux = np.loadtxt(heat_flux_file)
    print("  Interpolate bheatflx")
    new_bheatflx = scipy.interpolate.griddata(
        list(zip(heat_flux[:, 1], heat_flux[:, 0])),
        heat_flux[:, 2],
        (new.y_grid, new.x_grid),
        method="nearest",
    )
    nc_new.variables["bheatflx"][0, :, :] = new_bheatflx[:, :]
    nc_new.variables["bheatflx"].long_name = "basal heat flux"
    nc_new.variables[
        "bheatflx"
    ].standard_name = "upward_geothermal_heat_flux_at_ground_level_in_land_ice"
    nc_new.variables["bheatflx"].units = "mW m-2"
    nc_new.variables["bheatflx"].ancillary_variables = "bheatflxerr"
    nc_new.variables["bheatflx"].grid_mapping = "epsg_3031"
    nc_new.variables["bheatflx"].coordinates = "lon lat"
    nc_new.variables["bheatflx"].reference = (
        "Martos, Yasmina M (2017): Antarctic geothermal heat flux distribution "
        "and estimated Curie Depths, links to gridded files. PANGAEA, "
        "https://doi.org/10.1594/PANGAEA.882503."
    )
    nc_new.variables["bheatflx"].comments = (
        "Resampled from 15km grid using 2D nearest neighbor interpolation; "
        "polar sterographic projection true scaled latitude not specified in "
        "dataset -- asumed -71 deg. (EPSG:3031)"
    )
    # FIXME nc_new.variables['bheatflx'].missing_value = FIXME

    print("\nHeat Flux Uncertainty")
    heat_flux_err = np.loadtxt(heat_flux_unc_file)
    print("  Interpolate bheatflxerr")
    new_bheatflxerr = scipy.interpolate.griddata(
        list(zip(heat_flux_err[:, 1], heat_flux_err[:, 0])),
        heat_flux_err[:, 2],
        (new.y_grid, new.x_grid),
        method="nearest",
    )

    nc_new.createVariable("bheatflxerr", "f4", ("time", "y1", "x1",))
    nc_new.variables["bheatflxerr"][0, :, :] = new_bheatflxerr[:, :]
    nc_new.variables["bheatflxerr"].long_name = "basal heat flux uncertainty"
    nc_new.variables[
        "bheatflxerr"
    ].standard_name = "upward_geothermal_heat_flux_at_ground_level_in_land_ice standard_error"

    nc_new.variables["bheatflxerr"].units = "mW m-2"
    nc_new.variables["bheatflxerr"].grid_mapping = "epsg_3031"
    nc_new.variables["bheatflxerr"].coordinates = "lon lat"
    nc_new.variables["bheatflxerr"].reference = (
        "Martos, Yasmina M (2017): Antarctic geothermal heat flux distribution "
        "and estimated Curie Depths, links to gridded files. PANGAEA, "
        "https://doi.org/10.1594/PANGAEA.882503. "
    )

    nc_new.variables["bheatflxerr"].comments = (
        "Resampled from 15km grid using 2D nearest neighbor interpolation"
        "; polar sterographic projection true scaled latitude not "
        "specified in dataset -- asumed -71 deg. (EPSG:3031)"
    )

    # FIXME nc_new.variables['bheatflxerr'].missing_value = FIXME
    # NOTE: bheatflxerr has no missing values because of the NNI with griddata.
    # Using another method would likely change this
    print("\nWrite metadata")
    # -------------------------
    # actually provide metadata
    # -------------------------
    nc_new.title = "CISM-style input dataset for ice-sheet models"
    nc_new.history = "Created {} by J. Kennedy & M. Kelleher.".format(
        datetime.now().strftime("%c")
    )

    nc_new.institution = "Oak Ridge National Laboratory"
    nc_new.references = "See https://github.com/mkstratos/cism-data"
    nc_new.Conventions = "CF-1.6"

    nc_new.variables["time"].long_name = "time"
    nc_new.variables["time"].units = "common_years since 2008-01-01 00:00:00"
    nc_new.variables["time"].calendar = "365_day"
    nc_new.variables["time"].comments = (
        "The initial time here is an estimate of the nominal date for "
        "Rignot (2011) InSAR velocity data. Because this is a synthesis "
        "of datasets across many time periods, the inital date is "
        "inherently fuzzy and should be changed to suit your purposes."
    )

    nc_new.variables["y1"].long_name = "y-coordinate of projection"
    nc_new.variables["y1"].standard_name = "projection_y_coordinate"
    # nc_new.variables['y1'].axis = 'Y'
    nc_new.variables["y1"].units = "m"

    nc_new.variables["x1"].long_name = "x-coordinate of projection"
    nc_new.variables["x1"].standard_name = "projection_x_coordinate"
    # nc_new.variables['x1'].axis = 'X'
    nc_new.variables["x1"].units = "m"

    nc_new.createVariable("epsg_3031", "c")
    nc_new.variables["epsg_3031"].grid_mapping_name = "polar_stereographic"
    nc_new.variables["epsg_3031"].straight_vertical_longitude_from_pole = 0.0
    nc_new.variables["epsg_3031"].latitude_of_projection_origin = -90.0
    nc_new.variables["epsg_3031"].standard_parallel = -71.0
    nc_new.variables["epsg_3031"].false_easting = 0.0
    nc_new.variables["epsg_3031"].false_northing = 0.0
    nc_new.variables["epsg_3031"].scale_factor_at_projection_origin = 1.0
    nc_new.variables["epsg_3031"].reference_ellipsoid_name = "WGS84"
    nc_new.variables["epsg_3031"].horizontal_datum_name = "WGS84"
    nc_new.variables["epsg_3031"].units = "m"
    nc_new.variables["epsg_3031"].proj4_string = (
        "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_"
        "0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
    )

    print("\nGet Projections")
    epsg3412, epsg3031, _ = projections.antarctica()
    nc_new, new = projections.grid_center_latlons(
        nc_new, new, epsg3031, "epsg_3031", ("y1", "x1")
    )

    nc_new.variables[
        "acab"
    ].long_name = "water equivalent surface mass balance"
    nc_new.variables[
        "acab"
    ].standard_name = "land_ice_lwe_surface_specific_mass_balance"
    nc_new.variables["acab"].units = "mm year-1"
    nc_new.variables["acab"].grid_mapping = "epsg_3031"
    nc_new.variables["acab"].coordinates = "lon lat"
    nc_new.variables["acab"].source = "J. T. M. Lenaerts"
    nc_new.variables["acab"].reference = (
        "Lenaerts, J. T. M., M. R. vanden Broeke, W. J. van deBerg, "
        "E. vanMeijgaard, and P. Kuipers Munneke (2012), A new, high‐resolution"
        " surface mass balance map of Antarctica (1979–2010) based on regional "
        "atmospheric climate modeling, Geophys. Res. Lett., 39, L04501, "
        "doi:10.1029/2011GL050713."
    )
    nc_new.variables[
        "acab"
    ].comments = (
        "Mean 1979--2010 SMB; 2D linear interpolation of 27 km dataset."
    )
    # FIXME nc_new.variables['acab'].missing_value = FIXME

    nc_new.variables[
        "artm"
    ].long_name = "annual mean air temperature (2 meter)"
    nc_new.variables["artm"].standard_name = "air_temperature"
    nc_new.variables["artm"].units = "degree_Celsius"
    nc_new.variables["artm"].grid_mapping = "epsg_3031"
    nc_new.variables["artm"].coordinates = "lon lat"
    nc_new.variables["artm"].source = "J. T. M. Lenaerts"
    nc_new.variables["artm"].reference = (
        "Lenaerts, J. T. M., M. R. vanden Broeke, W. J. van deBerg, "
        "E. vanMeijgaard, and P. Kuipers Munneke (2012), A new, high‐resolution"
        " surface mass balance map of Antarctica (1979–2010) based on regional "
        "atmospheric climate modeling, Geophys. Res. Lett., 39, L04501, "
        "doi:10.1029/2011GL050713."
    )
    nc_new.variables[
        "artm"
    ].comments = (
        "Mean 1979--2010 t2m; 2D linear interpolation of 27 km dataset."
    )
    # FIXME nc_new.variables['artm'].missing_value = FIXME

    nc_new.variables["topg"].long_name = "bed topography"
    nc_new.variables["topg"].standard_name = "bedrock_altitude"
    nc_new.variables["topg"].units = "m"
    nc_new.variables["topg"].grid_mapping = "epsg_3031"
    nc_new.variables["topg"].coordinates = "lon lat"
    nc_new.variables["topg"].reference = (
        "Fretwell, P., et al.: Bedmap2: improved ice bed, surface and "
        "thickness datasets for Antarctica, The Cryosphere, 7, 375-393, "
        "https://doi.org/10.5194/tc-7-375-2013, 2013."
    )
    nc_new.variables[
        "topg"
    ].comments = (
        "Resampled from 5km grid using 2D nearest neighbor interpolation"
    )
    # FIXME nc_new.variables['topg'].missing_value = FIXME

    nc_new.variables["vx"].long_name = "surface x velocity"
    nc_new.variables["vx"].standard_name = "land_ice_surface_x_velocity"
    nc_new.variables["vx"].units = "m year-1"
    nc_new.variables["vx"].grid_mapping = "epsg_3031"
    nc_new.variables["vx"].coordinates = "lon lat"
    nc_new.variables["vx"].source = "NSIDC; Dataset ID: NSIDC-0484 v1.1"
    nc_new.variables["vx"].reference = (
        "Rignot, E., J. Mouginot, and B. Scheuchl. 2011. Ice Flow of the "
        "Antarctic Ice Sheet, Science. 333. 1427-1430. "
        "https://doi.org/10.1126/science.1208336. Mouginot J., B. Scheuchl and "
        "E. Rignot (2012), Mapping of Ice Motion in Antarctica Using "
        "Synthetic-Aperture Radar Data, Remote Sensing, doi 10.3390/rs4092753"
    )

    nc_new.variables[
        "vx"
    ].comments = "2D linear interpolation of provided 450 m dataset"
    # FIXME nc_new.variables['vx'].missing_value = FIXME

    nc_new.variables["vy"].long_name = "surface y velocity"
    nc_new.variables["vy"].standard_name = "land_ice_surface_y_velocity"
    nc_new.variables["vy"].units = "m year-1"
    nc_new.variables["vy"].grid_mapping = "epsg_3031"
    nc_new.variables["vy"].coordinates = "lon lat"
    nc_new.variables["vx"].source = "NSIDC; Dataset ID: NSIDC-0484 v1.1"
    nc_new.variables["vx"].reference = (
        "Rignot, E., J. Mouginot, and B. Scheuchl. 2011. Ice Flow of the "
        "Antarctic Ice Sheet, Science. 333. 1427-1430. "
        "https://doi.org/10.1126/science.1208336. Mouginot J., B. Scheuchl and "
        "E. Rignot (2012), Mapping of Ice Motion in Antarctica Using "
        "Synthetic-Aperture Radar Data, Remote Sensing, doi 10.3390/rs4092753"
    )
    nc_new.variables[
        "vx"
    ].comments = "2D linear interpolation of provided 450 m dataset"
    # FIXME nc_new.variables['vy'].missing_value = FIXME

    # FIXME This should really be vxerr and vyerr to follow CF conventions
    #   can use standard name modifier and associate vx and xv to them via an
    #   ancillary_variable attribute
    verr = np.ma.masked_values(nc_new.variables["verr"][0, :, :], 0.0)
    if FILL:
        nc_new.variables["verr"][0, :, :] = verr[:, :].filled(20.0)
    else:
        nc_new.variables["verr"][0, :, :] = verr[:, :]
    nc_new.variables[
        "verr"
    ].long_name = "magnitude of surface velocity error estimate"
    nc_new.variables["verr"].units = "m year-1"
    nc_new.variables["verr"].grid_mapping = "epsg_3031"
    nc_new.variables["verr"].coordinates = "lon lat"
    nc_new.variables["verr"].source = "NSIDC; Dataset ID: NSIDC-0484 v1.1"
    nc_new.variables["verr"].reference = (
        "Rignot, E., J. Mouginot, and B. Scheuchl. 2011. Ice Flow of "
        "the Antarctic Ice Sheet, Science. 333. 1427-1430. "
        "https://doi.org/10.1126/science.1208336.Mouginot J., "
        "B. Scheuchl and E. Rignot (2012), Mapping of Ice Motion in "
        "Antarctica Using Synthetic-Aperture Radar Data, Remote "
        "Sensing, doi 10.3390/rs4092753 "
    )
    nc_new.variables[
        "verr"
    ].comments = "2D linear interpolation of provided 450 m dataset"
    if not FILL:
        nc_new.variables["verr"].missing_value = -9999.0

    print("\nCryosat")
    nc_cryosat = ncfunc.get_nc_file(cryosat_file, "r")
    cryosat = projections.DataGrid()
    # NOTE: Grid coordinates relate to lower-left
    # of cell, not grid center like CISM
    cryosat.yll = nc_cryosat.variables["y1"][:]  # m
    cryosat.xll = nc_cryosat.variables["x1"][:]  # m
    cryosat.dy = cryosat.yll[1] - cryosat.yll[0]
    cryosat.dx = cryosat.xll[1] - cryosat.xll[0]
    cryosat.y = cryosat.yll + cryosat.dy / 2.0
    cryosat.x = cryosat.xll + cryosat.dx / 2.0
    cryosat.ny = cryosat.y.shape[0]
    cryosat.nx = cryosat.x.shape[0]
    cryosat.make_grid()

    cryosat.dhdt = np.ma.masked_values(
        nc_cryosat.variables["dzdt"][-1, :, :], -9999.0
    )
    cryosat.ygma = cryosat.y_grid[~cryosat.dhdt.mask]
    cryosat.xgma = cryosat.x_grid[~cryosat.dhdt.mask]
    cryosat.zgma = cryosat.dhdt[~cryosat.dhdt.mask]
    print("  Transform")
    cryosat.txgma, cryosat.tygma = pyproj.transform(
        epsg3412, epsg3031, cryosat.xgma.flatten(), cryosat.ygma.flatten()
    )
    print("  Interpolate")
    dhdt_interp = scipy.interpolate.griddata(
        list(zip(cryosat.tygma, cryosat.txgma)),
        cryosat.zgma.flatten(),
        (new.y_grid, new.x_grid),
        method="linear",
    )

    print("  Mask")
    msk_ice = nc_new.variables["thk"][-1, :, :] > np.nextafter(0.0, 1.0)
    msk_grounded = np.logical_and(
        msk_ice,
        (-1023.0 / 917.0) * nc_new.variables["topg"][-1, :, :]
        < np.nextafter(
            nc_new.variables["thk"][-1, :, :],
            nc_new.variables["thk"][-1, :, :] + 1.0,
        ),
    )
    msk_polehole = nc_new.variables["lat"] < np.nextafter(-88.0, 1.0)
    msk = np.logical_and(msk_grounded, msk_ice)
    msk = np.logical_and(msk, msk_polehole)

    dhdt = np.ma.masked_where(~msk_grounded, dhdt_interp)

    nc_new.createVariable("dhdt", "f4", ("time", "y1", "x1",))
    if FILL:
        nc_new.variables["dhdt"][0, :, :] = dhdt[:, :].filled(0.0)
    else:
        nc_new.variables["dhdt"][0, :, :] = dhdt[:, :]
    nc_new.variables["dhdt"].long_name = "observed thickness tendency"
    nc_new.variables["dhdt"].standard_name = "tendency_of_land_ice_thickness"
    nc_new.variables["dhdt"].units = "m year-1"
    nc_new.variables["dhdt"].ancillary_variables = "dhdterr"
    nc_new.variables["dhdt"].grid_mapping = "epsg_3031"
    nc_new.variables["dhdt"].coordinates = "lon lat"
    nc_new.variables["dhdt"].source = "M. McMillan and A. Shepherd"
    nc_new.variables["dhdt"].reference = (
        "Mcmillan, M., A. Shepherd, A. Sundal, K. Briggs, A. Muir, "
        "A. Ridout, A. Hogg, and D. Wingham, 2014: Increased ice losses "
        "from Antarctica detected by CryoSat-2. Geophys. Res. Lett, "
        "doi:10.1002/2014GL060111."
    )
    nc_new.variables["dhdt"].comments = (
        "As per the request of the authors (M. McMillan & A. Shepherd), "
        "these data are not to be shared outside of this project "
        "(PISCEES). They are to be used for optimization and model "
        "validation purposes only, as the original authors still have "
        "plans to use them for other studies of their own. They are ok "
        "with us using them for optimization and validation with the "
        "caveat that we should communicate further with them about their "
        "use prior to publishing any stuides that use them. Also, if the "
        "data are used for presentations, we should acknowldege the "
        "authors as the source of the data. For any further questions, "
        "please check with S. Price or D. Martin."
    )
    if not FILL:
        nc_new.variables["dhdt"].missing_value = -9999.0

    dhdterr = np.zeros(new.x_grid.shape)
    dhdterr[~msk] = 1.10  # m/year
    dhdterr[~msk] = 1.10  # m/year

    nc_new.createVariable("dhdterr", "f4", ("time", "y1", "x1",))
    nc_new.variables[
        "dhdt"
    ].long_name = "observed thickness tendency uncertainty"
    nc_new.variables[
        "dhdt"
    ].standard_name = "tendency_of_land_ice_thickness standard_error"
    nc_new.variables["dhdt"].units = "m year-1"
    nc_new.variables["dhdt"].grid_mapping = "epsg_3031"
    nc_new.variables["dhdt"].coordinates = "lon lat"
    nc_new.variables["dhdt"].source = "M. McMillan and A. Shepherd"
    nc_new.variables["dhdt"].reference = (
        "Mcmillan, M., A. Shepherd, A. Sundal, K. Briggs, A. Muir, "
        "A. Ridout, A. Hogg, and D. Wingham, 2014: Increased ice "
        "losses from Antarctica detected by CryoSat-2. "
        "Geophys. Res. Lett, doi:10.1002/2014GL060111."
    )
    nc_new.variables["dhdt"].comments = (
        "Uncertainty estimates over the ice sheet proper, exlcuding the "
        "pole-hole and ice shelves, represent the std deviation of the "
        "difference between CryoSat and airborn altimetry. Over the pole-hole "
        "and ice-shelves where data is missing, we've doubled the uncertainty."
    )

    nc_new.close()
    nc_berr.close()
    nc_bike.close()
    nc_subm.close()

    if coarse_out:
        finalize.coarsen(args, "epsg_3031", f_1km_out, f_template, coarse_list)
    print("--" * 5 + "DONE" + "--" * 5)


if __name__ == "__main__":
    main(coarse_out=False)
