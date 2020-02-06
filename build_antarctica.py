#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build CISM style input file from multiple sources using xarray.
"""
import os
from datetime import datetime
import numpy as np
import scipy
import scipy.interpolate
import xarray as xr
from nco import Nco
import pyproj
from util import projections

__author__ = "Michael Kelleher"


def xr_load(in_file, **kwargs):
    """Load netCDF data."""
    return xr.open_dataset(in_file)


def np_load(in_file, **kwargs):
    """Load text file."""
    return np.loadtxt(in_file)


def load_cryosat(in_file, **kwargs):
    """Load and process Cryosat2 data."""
    grid_from = kwargs["grid_from"]
    grid_to = kwargs["grid_to"]
    cfg = kwargs["cfg"]
    in_var = kwargs["in_var"]
    # NOTE: Grid coordinates relate to lower-left of cell,
    # not grid center like CISM
    din = xr.open_dataset(in_file)
    _xin = din[cfg["coords"]["x"]]
    _yin = din[cfg["coords"]["y"]]
    x_in_grid, y_in_grid = np.meshgrid(_xin, _yin, indexing="ij")
    mask_data = np.ma.masked_values(din[in_var], -9999.0)
    x_in_grid = x_in_grid[~mask_data.mask]
    y_in_grid = y_in_grid[~mask_data.mask]

    x_transform, y_transform = pyproj.transform(
        grid_from, grid_to, x_in_grid, y_in_grid
    )
    return np.vstack(x_transform, y_transform, mask_data.flatten()).T


def load_hf(in_file, **kwargs):
    """"Load Basal Heat Flux data, create gridded xarray.DataArray."""
    assert "cfg" in kwargs, "Missing required argument cfg"
    assert "in_var" in kwargs, "Missing required argument in_var"
    cfg = kwargs["cfg"]
    in_var = kwargs["in_var"]
    d_x = cfg.get("dx", 15_000)
    d_y = cfg.get("dy", 15_000)
    data = np_load(in_file)
    xin = data[:, 0]
    yin = data[:, 1]
    hf_in = data[:, 2]
    x_ext = np.arange(xin.min(), xin.max() + d_x, d_x)
    y_ext = np.arange(yin.min(), yin.max() + d_y, d_y)
    hf_grid = np.zeros((y_ext.shape[0], x_ext.shape[0]))

    for idx in range(hf_in.shape[0]):
        ii = np.where(x_ext == xin[idx])[0]
        jj = np.where(y_ext == yin[idx])[0]
        hf_grid[jj, ii] = hf_in[idx]

    hf_grid = xr.DataArray(
        hf_grid,
        coords={cfg["coords"]["x"]: x_ext, cfg["coords"]["y"]: y_ext},
        dims=[cfg["coords"]["y"], cfg["coords"]["x"]],
    )
    return xr.Dataset({in_var: hf_grid})


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
    x_grid, y_grid = np.meshgrid(xout, yout, indexing="ij")
    new_data = np.zeros(x_grid.shape)

    for _ii in range(0, y_grid.shape[1]):
        new_data[:, _ii] = to_new.ev(y_grid[:, _ii], x_grid[:, _ii])

    return new_data


def interp(in_cfg, in_var, out_cfg, output):
    """Take input data, interpolate, assgin metadata, return DataArray."""
    print(f"Interpolating {in_var} from {in_cfg['file']}")
    in_data = in_cfg["load"](in_cfg["file"], cfg=in_cfg, in_var=in_var)
    xout = output[out_cfg["coords"]["x"]]
    yout = output[out_cfg["coords"]["y"]]

    if isinstance(in_data, xr.Dataset):
        nx_in = in_data[in_cfg["coords"]["x"]].shape[0]
        ny_in = in_data[in_cfg["coords"]["y"]].shape[0]
        if nx_in >= xout.shape[0] and ny_in >= yout.shape[0]:
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
        out_grid = np.meshgrid(xout, yout, indexing="ij")
        intp_data = scipy.interpolate.griddata(
            in_data[:, :2][:, ::-1],
            in_data[:, 2],
            tuple(out_grid),
            method="nearest",
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
        # Assign variable metadata
        _outdata = intp_data.assign_attrs(in_cfg["meta"][in_var])

        # Assign common dataset metadata, if it exists
        if in_cfg.get("cmeta"):
            _outdata = _outdata.assign_attrs(in_cfg["cmeta"])

    return _outdata.expand_dims("time", 0)


def output_setup(input_file, output_file):
    """Create output file by copying old 1 km grid file."""
    nco = Nco()
    nco.ncks(
        input=input_file,
        output=output_file,
        options=["-x", "-v", "acab_alb,artm_alb,dzdt"],
    )
    return xr.open_dataset(output_file).load()


def grid_center_lat_lon(dset, proj, proj_var_name, cvars=("y", "x")):
    """Create latitude/longitude variables based on proj4 projection."""
    lon_attrs = {
        "long_name": "grid center longitude",
        "standard_name": "longitude",
        "units": "degrees_east",
        "grid_mapping": proj_var_name,
    }
    lat_attrs = {
        "long_name": "grid center latitude",
        "standard_name": "latitude",
        "units": "degrees_north",
        "grid_mapping": proj_var_name,
    }
    x_grid, y_grid = np.meshgrid(dset[cvars[1]], dset[cvars[0]])
    lon_grid, lat_grid = proj(x_grid.ravel(), y_grid.ravel(), inverse=True)

    lon_grid = lon_grid.reshape(x_grid.shape)
    lat_grid = lat_grid.reshape(y_grid.shape)

    dset["lon"] = xr.DataArray(
        lon_grid,
        dims=cvars,
        coords={cvar: dset[cvar] for cvar in cvars},
        attrs=lon_attrs,
    )
    dset["lat"] = xr.DataArray(
        lat_grid,
        dims=cvars,
        coords={cvar: dset[cvar] for cvar in cvars},
        attrs=lat_attrs,
    )


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
                },
                "melt_steadystate": {
                    "long_name": "steady state sub-shelf melt rate",
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
            },
            "cmeta": {
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
                },
                "errbed": {
                    "long_name": "ice thickness error",
                    "standard_name": "land_ice_thickness standard_error",
                    "units": "m",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
            },
            "cmeta": {
                "source": "BedMachine Antarctica, Mathieu Morlighem",
                "reference": (
                    "Morlighem M. et al., (2019), Deep glacial troughs and "
                    "stabilizing ridges unveiled beneath the margins of the "
                    "Antarctic ice sheet, Nature Geoscience (accepted), "
                    "doi:10.1038/s41561-019-0510-8"
                ),
                "comments": (
                    "Obtained from NSIDC: https://nsidc.org/data/nsidc-0756. "
                    "Resampled from 500m grid using Nearest Neighbor -> 1km"
                ),
            },
        },
        "heatflux": {
            "file": "data/Martos-AIS-heatFlux/Antarctic_GHF.xyz",
            "load": load_hf,
            "vars": ["bheatflx"],
            "coords": {"x": "x", "y": "y"},
            "dx": 15000,
            "dy": 15000,
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
                }
            },
            "cmeta": {
                "reference": (
                    "Martos, Yasmina M (2017): Antarctic geothermal heat flux "
                    "distribution and estimated Curie Depths, links to gridded "
                    "files. PANGAEA, https://doi.org/10.1594/PANGAEA.882503."
                ),
                "comments": (
                    "Resampled from 15km grid using 2D nearest neighbor "
                    "interpolation; polar sterographic projection true scaled "
                    "latitude not specified in dataset -- assumed 71 deg. "
                    "(EPSG:3031)"
                ),
            },
        },
        "heatflux_unc": {
            "file": "data/Martos-AIS-heatFlux/Antarctic_GHF_uncertainty.xyz",
            "load": load_hf,
            "vars": ["bheatflxerr"],
            "coords": {"x": "x", "y": "y"},
            "dx": 15000,
            "dy": 15000,
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
                }
            },
        },
        "cryosat": {
            "file": "data/Cryosat2/CS2_dzdt.nc",
            "load": xr_load,
            "vars": ["dzdt", "dzdterr"],
            "coords": {"x": "x1", "y": "y1"},
            "meta": {
                "dzdt": {
                    "long_name": "observed thickness tendency",
                    "standard_name": "tendency_of_land_ice_thickness",
                    "units": "m year-1",
                    "ancillary_variables": "dhdterr",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
                "dzdterr": {
                    "long_name": "observed thickness tendency uncertainty",
                    "standard_name": (
                        "tendency_of_land_ice_thickness standard error"
                    ),
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
            },
            "cmeta": {
                "source": "M. McMillan and A. Shepherd",
                "reference": (
                    "Mcmillan, M., A. Shepherd, A. Sundal, K. Briggs, A. Muir, "
                    "A. Ridout, A. Hogg, and D. Wingham, 2014: Increased ice "
                    "losses from Antarctica detected by CryoSat-2. Geophys. "
                    "Res. Lett, doi:10.1002/2014GL060111."
                ),
                "comments": (
                    "As per the request of the authors (M. McMillan & "
                    "A. Shepherd), these data are not to be shared outside of "
                    "this project (PISCEES). They are to be used for "
                    "optimization and model validation purposes only, as the "
                    "original authors still have plans to use them for other "
                    "studies of their own. They are ok with us using them for "
                    "optimization and validation with the caveat that we "
                    "should communicate further with them about their use "
                    "prior to publishing any stuides that use them. Also, if "
                    "the data are used for presentations, we should "
                    "acknowldege the authors as the source of the data. "
                    "For any further questions, please check with "
                    "S. Price or D. Martin."
                ),
            },
        },
        "smb": {
            "vars": ["acab", "artm"],
            "meta": {
                "acab": {
                    "long_name": "water equivalent surface mass balance",
                    "standard_name": (
                        "land_ice_lwe_surface_specific_mass_balance"
                    ),
                    "units": "mm year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
                "artm": {
                    "long_name": "annual mean air temperature (2 meter)",
                    "standard_name": "air_temperature",
                    "units": "degree_Celsius",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
            },
            "cmeta": {
                "source": "J. T. M. Lenaerts",
                "reference": (
                    "Lenaerts, J. T. M., M. R. vanden Broeke, W. J. van deBerg,"
                    " E. vanMeijgaard, and P. Kuipers Munneke (2012), A new, "
                    "high‐resolution surface mass balance map of Antarctica "
                    "(1979–2010) based on regional atmospheric climate "
                    "modeling, Geophys. Res. Lett., 39, L04501, "
                    "doi:10.1029/2011GL050713."
                ),
                "comments": (
                    "Mean 1979--2010 t2m; 2D linear interpolation of 27 km "
                    "dataset."
                ),
            },
        },
        "topg": {
            "vars": ["topg"],
            "meta": {
                "topg": {
                    "long_name": "bed topography",
                    "standard_name": "bedrock_altitude",
                    "units": "m",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
            },
            "cmeta": {
                "reference": (
                    "Fretwell, P., et al.: Bedmap2: improved ice bed, surface "
                    "and thickness datasets for Antarctica, The Cryosphere, 7, "
                    "375-393, https://doi.org/10.5194/tc-7-375-2013, 2013."
                ),
                "comments": (
                    "Resampled from 5km grid using 2D nearest "
                    "neighbor interpolation"
                ),
            },
        },
        "veloc": {
            "vars": ["vx", "vy", "verr"],
            "meta": {
                "vx": {
                    "long_name": "surface x velocity",
                    "standard_name": "land_ice_surface_x_velocity",
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "source": "NSIDC; Dataset ID: NSIDC-0484 v1.1",
                    "reference": (
                        "Rignot, E., J. Mouginot, and B. Scheuchl. 2011. "
                        "Ice Flow of the Antarctic Ice Sheet, Science. 333. "
                        "1427-1430. https://doi.org/10.1126/science.1208336. "
                        "Mouginot J., B. Scheuchl and E. Rignot (2012), "
                        "Mapping of Ice Motion in Antarctica Using "
                        "Synthetic-Aperture Radar Data, Remote Sensing, "
                        "doi 10.3390/rs4092753"
                    ),
                    "comments": (
                        "2D linear interpolation of provided 450 m dataset"
                    ),
                },
                "vy": {
                    "long_name": "surface y velocity",
                    "standard_name": "land_ice_surface_y_velocity",
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                    "source": "NSIDC; Dataset ID: NSIDC-0484 v1.1",
                    "reference": (
                        "Rignot, E., J. Mouginot, and B. Scheuchl. 2011. "
                        "Ice Flow of the Antarctic Ice Sheet, Science. 333. "
                        "1427-1430. https://doi.org/10.1126/science.1208336. "
                        "Mouginot J., B. Scheuchl and E. Rignot (2012), "
                        "Mapping of Ice Motion in Antarctica Using "
                        "Synthetic-Aperture Radar Data, Remote Sensing, "
                        "doi 10.3390/rs4092753"
                    ),
                    "comments": (
                        "2D linear interpolation of provided 450 m dataset"
                    ),
                },
                "verr": {
                    "long_name": "magnitude of surface velocity error estimate",
                    "units": "m year-1",
                    "grid_mapping": "epsg_3031",
                    "coordinates": "lon lat",
                },
            },
            "cmeta": {
                "source": "NSIDC; Dataset ID: NSIDC-0484 v1.1",
                "reference": (
                    "Rignot, E., J. Mouginot, and B. Scheuchl. 2011. Ice Flow "
                    "of the Antarctic Ice Sheet, Science. 333. 1427-1430. "
                    "https://doi.org/10.1126/science.1208336. Mouginot J., "
                    "B. Scheuchl and E. Rignot (2012), Mapping of Ice Motion "
                    "in Antarctica Using Synthetic-Aperture Radar Data, Remote "
                    "Sensing, doi 10.3390/rs4092753"
                ),
                "comments": (
                    "2D linear interpolation of provided 450 m dataset"
                ),
            },
        },
    }
    # Get metadata from heatflux config (so it doesn't have to be in here twice)
    input_config["heatflux_unc"]["cmeta"] = input_config["heatflux"]["cmeta"]

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
                "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0"
                " +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
            ),
        },
    }

    output_file = (
        f"ncs/antarctica_1km_{datetime.now().strftime('%Y_%m_%d')}.nc"
    )

    output_cfg = {"coords": {"x": "x1", "y": "y1"}}
    output = output_setup(input_config["1km_in"]["file"], output_file)

    # Unlink from file, since output_setup loads netCDF to memory
    # this allows writing back to the same file we open
    output.close()

    output = output.assign_attrs(output_metadata)
    for variable in output_variables:
        if variable == "epsg_3031":
            output[variable] = b""
            output[variable] = output[variable].assign_attrs()

        if variable in output:
            output[variable] = output[variable].assign_attrs(
                output_variables[variable]
            )

    # Map required output variable to a dataset and input variable
    # (OUTPUT VARIABLE, INPUT DATASET, INPUT VARIABLE NAME)
    inout_map = [
        ("thk", "bedmachine", "thickness"),
        ("thkerr", "bedmachine", "errbed"),
        ("bheatflx", "heatflux", "bheatflux"),
        ("bheatflxerr", "heatflux_unc", "bheatflxerr"),
        ("subm", "rignot_subshelf", "melt_actual"),
        ("subm_ss", "rignot_subshelf", "melt_steadystate"),
    ]

    all_exist = True

    # Add metadata for extant variables
    ext_vars = [
        ("smb", "acab"),
        ("smb", "artm"),
        ("topg", "topg"),
        ("veloc", "verr"),
        ("veloc", "vx"),
        ("veloc", "vy"),
    ]
    for dset, var in ext_vars:
        output[var] = output[var].assign_attrs(input_config[dset]["meta"][var])
        output[var] = output[var].assign_attrs(input_config[dset]["cmeta"])

    for _, ds_in, in_var in inout_map:
        if not os.path.exists(input_config[ds_in]["file"]):
            all_exist = False
            print(f"MISSING FILE: {input_config[ds_in]['file']}")
        else:
            print(f"FILE FOUND: {input_config[ds_in]['file']}")

    epsg3412, epsg3031, _ = projections.antarctica()
    grid_center_lat_lon(output, epsg3031, "epsg_3031", cvars=("y1", "x1"))

    print("-" * 50)
    if all_exist:
        for out_var, ds_in, in_var in inout_map:
            output[out_var] = interp(
                input_config[ds_in], in_var, output_cfg, output
            )

    # TODO: Check masking
    output.to_netcdf(output_file)


if __name__ == "__main__":
    main()
