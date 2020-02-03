#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build CISM style input file from multiple sources using xarray.
"""
from datetime import datetime
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import xarray as xr
from nco import Nco

__author__ = "Michael Kelleher"


def xr_load(in_file):
    """Load netCDF data."""
    return xr.open_dataset(in_file)


def np_load(in_file):
    """Load text file."""
    return np.loadtxt(in_file)


def rect_bivariate_interp(xin, yin, data, xout, yout):
    """Interpolate gridded data on rectangular grid with Bivariate spline."""
    if yin[1] - yin[0] < 0:
        _yslice = slice(None, None, -1)
    else:
        _yslice = slice(None, None, None)

    if xin[1] - xin[0] < 0:
        _xslice = slice(None, None, -1)
    else:
        _xslice = slice(None, None, None)

    to_new = scipy.interpolate.RectBivariateSpline(
        yin[_yslice], xin[_xslice], data[_yslice, _xslice], kx=1, ky=1, s=0
    )
    x_grid, y_grid = np.meshgrid(xout, yout)
    new_data = np.zeros(x_grid.shape)

    for ii in range(0, y_grid.shape[1]):
        new_data[:, ii] = to_new.ev(y_grid[:, ii], x_grid[:, ii])

    return new_data


def interp(in_cfg, in_var, out_cfg, output):
    """Take input data, interpolate, assgin metadata, return DataArray."""
    print(f"Interpolating {in_var} from {in_cfg['file']}")
    in_data = in_cfg["load"](in_cfg["file"])
    xout = output[out_cfg["coords"]["x"]]
    yout = output[out_cfg["coords"]["y"]]
    nx_out = xout.shape[0]
    ny_out = yout.shape[0]

    if isinstance(in_data, xr.Dataset):
        nx_in = in_data[in_cfg["coords"]["x"]].shape[0]
        ny_in = in_data[in_cfg["coords"]["y"]].shape[0]
        if nx_in >= nx_out and ny_in >= ny_out:
            # Use nearest neighbour when input data is higher res than output
            print("    using nearest neighbour")
            intp_data = in_data[in_var].interp(
                **{in_cfg["coords"]["x"]: xout, in_cfg["coords"]["y"]: yout},
                method="nearest",
            )
            for coord in ["y", "x"]:
                if in_cfg["coords"][coord] != out_cfg["coords"][coord]:
                    # Drop original coordinate copy from the output DataArray
                    intp_data = intp_data.drop(in_cfg["coords"][coord])
        else:
            # Otherwise use Bivariate Spline on rectangular grid
            print("    using bivariate spline")
            xin = in_data[in_cfg["coords"]["x"]]
            yin = in_data[in_cfg["coords"]["y"]]
            intp_data = rect_bivariate_interp(
                xin.values,
                yin.values,
                in_data[in_var].values,
                xout.values,
                yout.values,
            )
    else:
        xin = in_data[:, 0]
        yin = in_data[:, 1]
        zin = in_data[:, 2]
        out_grid = np.meshgrid(xout, yout)
        intp_data = scipy.interpolate.griddata(
            list(zip(yin, xin)), zin, out_grid, method="nearest"
        )

    if not isinstance(intp_data, xr.DataArray):
        _outdata = xr.DataArray(
            intp_data,
            coords={
                out_cfg["coords"]["x"]: xout,
                out_cfg["coords"]["y"]: yout,
            },
            dims=[out_cfg["coords"]["y"], out_cfg["coords"]["x"]],
            attrs=in_cfg["meta"][in_var],
        )
    else:
        _outdata = intp_data.assign_attrs(in_cfg["meta"][in_var])

    return _outdata.expand_dims("time", 0)


def output_setup(input_file, output_file):
    """Create output file by copying old 1 km grid file."""
    nco = Nco()
    nco.ncks(
        input=input_file,
        output=output_file,
        options=["-x", "-v", "acab_alb,artm_alb,dzdt"],
    )
    return xr.open_dataset(output_file)


def main():
    """Define configuration and required variables."""
    input_config = {
        "1km_in": {
            "file": "ncs/antarctica_1km_2017_05_03.nc",
            "vars": ["acab_alb", "artm_alb", "dzdt"],
            "coords": {"x": "x1", "y": "y1"},
        },
        # Rignot Subshelf Melt rates file
        "rignot_subshelf": {
            "file": (
                "data/Rignot-subShelfMeltRates/"
                "Ant_MeltingRate.flipNY.newAxes.nc"
            ),
            "load": xr_load,
            "vars": ["melt_actual", "melt_steadystate"],
            "coords": {"x": "x1", "y": "y1"},
            "meta": {
                "melt_actual": {
                    "long_name": "sub-shelf melt rate",
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "reference": (
                        "Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, "
                        "2013: Ice-Shelf Melting Around Antarctica. Science, "
                        "341, 266-270, doi:10.1126/science.1235798."
                    ),
                    "source": "J. Mouginot",
                    "comments": (
                        "2D linear interpolation of provided dataset. Authors "
                        "request we communicate with them prior to publishing "
                        "and acknowledge them as the source of the data in "
                        "presentations."
                    ),
                },
                "melt_steadystate": {
                    "long_name": "steady state sub-shelf melt rate",
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "reference": (
                        "Rignot, E., S. Jacobs, J. Mouginot, and B. Scheuchl, "
                        "2013: Ice-Shelf Melting Around Antarctica. Science, "
                        "341, 266-270, doi:10.1126/science.1235798."
                    ),
                    "source": "J. Mouginot",
                    "comments": (
                        "2D linear interpolation of provided dataset. Authors "
                        "request we communicate with them prior to publishing "
                        "and acknowledge them as the source of the data in "
                        "presentations."
                    ),
                },
            },
        },
        # BedMachine thickness, topography
        "bedmachine": {
            "file": (
                "data/500m.MassConsBed.AIS.Morlighem.2019/"
                "BedMachineAntarctica_2019-11-05_v01.nc"
            ),
            "load": xr_load,
            "vars": ["thickness", "errbed"],
            "coords": {"x": "x", "y": "y"},
            "meta": {
                "thickness": {
                    "long_name": "ice thickness",
                    "standard_name": "land_ice_thickness",
                    "units": "m",
                    "ancillary_variables": "thkerr",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "source": "BedMachine Antarctica, Mathieu Morlighem",
                    "reference": (
                        "Morlighem M. et al., (2019), Deep glacial troughs "
                        "and stabilizing ridges unveiled beneath the margins "
                        "of the Antarctic ice sheet, Nature Geoscience ("
                        "accepted), doi:10.1038/s41561-019-0510-8"
                    ),
                    "comments": (
                        "Obtained from NSIDC: "
                        "https://nsidc.org/data/nsidc-0756. "
                        "Resampled from 500m grid using Nearest Neighbor -> 1km"
                    ),
                },
                "errbed": {
                    "long_name": "ice thickness error",
                    "standard_name": "land_ice_thickness standard_error",
                    "units": "m",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "source": "BedMachine Antarctica, Mathieu Morlighem",
                    "reference": (
                        "Morlighem M. et al., (2019), Deep glacial troughs "
                        "and stabilizing ridges unveiled beneath the margins "
                        "of the Antarctic ice sheet, Nature Geoscience ("
                        "accepted), doi:10.1038/s41561-019-0510-8"
                    ),
                    "comments": (
                        "Obtained from NSIDC: "
                        "https://nsidc.org/data/nsidc-0756. "
                        "Resampled from 500m grid using Nearest Neighbor -> 1km"
                    ),
                },
            },
        },
        "heatflux": {
            "file": "data/Martos-AIS-heatFlux/Antarctic_GHF.xyz",
            "load": np_load,
            "vars": ["bheatflx"],
            "coords": {"x": "x", "y": "y"},
            "meta": {
                "bheatflux": {
                    "long_name": "basal heat flux",
                    "standard_name": (
                        "upward_geothermal_heat_flux_at_ground_level_in_land_ice"
                    ),
                    "units": "mW m-2",
                    "ancillary_variables": "bheatflxerr",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "reference": (
                        "Martos, Yasmina M (2017): Antarctic geothermal heat "
                        "flux distribution and estimated Curie Depths, links "
                        "to gridded files. PANGAEA, "
                        "https://doi.org/10.1594/PANGAEA.882503."
                    ),
                    "comments": (
                        "Resampled from 15km grid using 2D nearest neighbor "
                        "interpolation; polar sterographic projection true "
                        "scaled latitude not specified in dataset -- assumed "
                        "-71 deg. (EPSG:3031)"
                    ),
                }
            },
        },
        "heatflux_unc": {
            "file": "data/Martos-AIS-heatFlux/Antarctic_GHF_uncertainty.xyz",
            "load": np_load,
            "vars": ["bheatflxerr"],
            "meta": {
                "bheatflxerr": {
                    "long_name": "basal heat flux uncertainty",
                    "standard_name": (
                        "upward_geothermal_heat_flux_at_ground_"
                        "level_in_land_ice standard_error"
                    ),
                    "units": "mW m-2",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "reference": (
                        "Martos, Yasmina M (2017): Antarctic geothermal heat "
                        "flux distribution and estimated Curie Depths, links "
                        "to gridded files. PANGAEA, "
                        "https://doi.org/10.1594/PANGAEA.882503."
                    ),
                    "comments": (
                        "Resampled from 15km grid using 2D nearest "
                        "neighbor interpolation; polar sterographic projection "
                        "true scaled latitude not specified in dataset -- "
                        "assumed 71 deg. (EPSG:3031)"
                    ),
                }
            },
        },
        "cryosat": {
            "file": "data/Cryosat2/CS2_dzdt.nc",
            "load": xr_load,
            "vars": ["dzdt"],
            "coords": {"x": "x1", "y": "y1"},
            "meta": {
                "dzdt": {
                    "long_name": "observed thickness tendency",
                    "standard_name": "tendency_of_land_ice_thickness",
                    "units": "m year-1",
                    "ancillary_variables": "dhdterr",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "source": "M. McMillan and A. Shepherd",
                    "reference": (
                        "Mcmillan, M., A. Shepherd, A. Sundal, K. Briggs, "
                        "A. Muir, A. Ridout, A. Hogg, and D. Wingham, 2014: "
                        "Increased ice losses from Antarctica detected by "
                        "CryoSat-2. Geophys. Res. Lett, "
                        "doi:10.1002/2014GL060111."
                    ),
                    "comments": (
                        "As per the request of the authors (M. McMillan & A. "
                        "Shepherd), these data are not to be shared outside "
                        "of this project (PISCEES). They are to be used for "
                        "optimization and model validation purposes only, as "
                        "the original authors still have plans to use them "
                        "for other studies of their own. They are ok with us "
                        "using them for optimization and validation with the "
                        "caveat that we should communicate further with them "
                        "about their use prior to publishing any stuides that "
                        "use them. Also, if the data are used for "
                        "presentations, we should acknowldege the authors as "
                        "the source of the data. For any further questions, "
                        "please check with S. Price or D. Martin."
                    ),
                }
            },
        },
    }

    output_metadata = {
        "title": "CISM-style input dataset for ice-sheet models",
        "history": "Created {} by J. Kennedy & M. Kelleher.".format(
            datetime.now().strftime("%c")
        ),
        "institution": "Oak Ridge National Laboratory",
        "references": "See https://github.com/mkstratos/cism-data",
        "Conventions": "CF-1.6",
    }
    output_variables = {
        "time": {
            "long_name": "time",
            "units": "common_years since 2008-01-01 00:00:00",
            "calendar": "365_day",
            "comments": (
                "The initial time here is an estimate of the nominal date for "
                "Rignot (2011) InSAR velocity data. Because this is a "
                "synthesis of datasets across many time periods, the inital "
                "date is inherently fuzzy and should be changed to suit your "
                "purposes."
            ),
        },
        "y1": {
            "long_name": "y-coordinate of projection",
            "standard_name": "projection_y_coordinate",
            "axis": "Y",
            "units": "m",
        },
        "x1": {
            "long_name": "x-coordinate of projection",
            "standard_name": "projection_x_coordinate",
            "axis": "X",
            "units": "m",
        },
        "epsg_3031": {
            "grid_mapping_name": "polar_stereographic",
            "straight_vertical_longitude_from_pole": 0.0,
            "latitude_of_projection_origin": -90.0,
            "standard_parallel": -71.0,
            "false_easting": 0.0,
            "false_northing": 0.0,
            "scale_factor_at_projection_origin": 1.0,
            "reference_ellipsoid_name": "WGS84",
            "horizontal_datum_name": "WGS84",
            "units": "m",
            "proj4_string": (
                "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_"
                "0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
            ),
        },
    }
    output_file = (
        f"ncs/antarctica_1km_{datetime.now().strftime('%Y_%m_%d')}.nc"
    )
    output_cfg = {"coords": {"x": "x1", "y": "y1"}}
    output = output_setup(input_config["1km_in"]["file"], output_file)
    output = output.assign_attrs(output_metadata)
    for variable in output_variables:
        if variable == "epsg_3031":
            output[variable] = b""
            output[variable] = output[variable].assign_attrs()

        if variable in output:
            output[variable] = output[variable].assign_attrs(
                output_variables[variable]
            )

    inout_map = [
        ("bedmachine", "thickness", "thk"),
        ("rignot_subshelf", "melt_steadystate", "subm_ss"),
    ]
    for ds_in, in_var, out_var in inout_map:
        output[out_var] = interp(
            input_config[ds_in], in_var, output_cfg, output
        )

    # TODO: Check masking
    # TODO: Setup epsg_3031 grid, lat lon coord variables
    # TODO: Write to netCDF


if __name__ == "__main__":
    main()
