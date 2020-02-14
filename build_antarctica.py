#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build CISM style input file from multiple sources using xarray.
"""
import importlib as ilib
import sys
from pathlib import Path
import json
from datetime import datetime
import numpy as np
import scipy
import scipy.interpolate
import xarray as xr
from nco import Nco
import pyproj
from util import projections


__author__ = "Michael Kelleher"

DATA_ROOT = Path("data")


def xr_load(in_file, **kwargs):
    """Load netCDF data."""
    return xr.open_dataset(in_file)


def np_load(in_file, **kwargs):
    """Load text file."""
    return np.loadtxt(in_file)


def load_cryosat(in_file, **kwargs):
    """Load and process Cryosat2 data."""
    cfg = kwargs["cfg"]
    grid_from = cfg["grid_from"]
    grid_to = cfg["grid_to"]
    in_var = kwargs["in_var"]
    # NOTE: Grid coordinates relate to lower-left of cell,
    # not grid center like CISM
    din = xr_load(in_file)
    _xin = din[cfg["coords"]["x"]]
    _yin = din[cfg["coords"]["y"]]
    x_in_grid, y_in_grid = np.meshgrid(_xin, _yin)
    x_transform, y_transform = pyproj.transform(
        grid_from, grid_to, x_in_grid, y_in_grid
    )
    x_transform = x_transform.mean(axis=0)
    y_transform = y_transform.mean(axis=1)

    data = din[in_var].values
    _out = xr.DataArray(
        np.ma.masked_values(data[-1, :, :], -9999),
        dims=[cfg["coords"]["y"], cfg["coords"]["x"]],
        coords={
            cfg["coords"]["y"]: y_transform,
            cfg["coords"]["x"]: x_transform,
        },
    )
    return xr.Dataset({in_var: _out})


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


def build_base_greenland(
    out_file, epsg_config, mapping_name, coords, d_meters=1000
):
    """Build base coordinates for EPSG:3413 projection."""

    proj_epsg3413, proj_mcb = projections.greenland()
    projs = {
        "epsg_3413": {
            "prj": proj_epsg3413,
            "lat_0": 90.0,
            "lon_0": -45.0,
            "lat_ts": 70.0,
            "scale": 1.0,
        },
        "mcb": {
            "prj": proj_mcb,
            "lat_0": 90.0,
            "lon_0": 321.0,
            "lat_ts": 71.0,
            "scale": 1.0,
        },
    }

    with open(epsg_config, "r") as f:
        cfg = json.load(f)

    x1d = np.arange(cfg["ll"][0], cfg["ur"][0], d_meters)
    y1d = np.arange(cfg["ll"][1], cfg["ur"][1], d_meters)

    grid_attrs = {
        "grid_mapping_name": "polar_stereographic",
        "latitude_of_projection_origin": projs[mapping_name]["lon_0"],
        "straight_vertical_longitude_from_pole": projs[mapping_name]["lat_0"],
        "standard_parallel": projs[mapping_name]["lat_ts"],
        "proj_scale_factor": projs[mapping_name]["scale"],
        "false_easting": 0.0,
        "false_northing": 0.0,
        "ellipsoid": "WGS84",
        "datum": "WGS84",
        "units": "meters",
        "proj4_string": projections.proj_string(projs[mapping_name]["prj"]),
    }

    attrs = {
        coords[axis]: {
            "long_name": f"{axis}-coordinate in projection",
            "standard_name": f"projection_{axis}_coordinate",
            "axis": axis.capitalize(),
            "units": "m",
            "grid_mapping": mapping_name,
        }
        for axis in coords
    }

    epsg_grid = xr.Dataset({coords["x"]: x1d, coords["y"]: y1d})
    grid_map = xr.Variable(dims=[], data=None, attrs=grid_attrs)

    for axis in coords:
        epsg_grid[coords[axis]] = epsg_grid[coords[axis]].assign_attrs(
            attrs[coords[axis]]
        )
    epsg_grid[mapping_name] = grid_map

    epsg_grid.to_netcdf(out_file)
    return epsg_grid


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
            method = "nearest"
        else:
            # Otherwise use linear
            method = "linear"

        print(f"    using {method}")
        intp_data = in_data[in_var].interp(
            **{in_cfg["coords"]["x"]: xout, in_cfg["coords"]["y"]: yout},
            method=method,
        )
        for coord in ["y", "x"]:
            if in_cfg["coords"][coord] != out_cfg["coords"][coord]:
                # Drop original coordinate copy from the output DataArray
                intp_data = intp_data.drop(in_cfg["coords"][coord])
        # else:
        #     # Otherwise use Bivariate Spline on rectangular grid
        #     print("    using bivariate spline")
        #     xin = in_data[in_cfg["coords"]["x"]]
        #     yin = in_data[in_cfg["coords"]["y"]]
        #     intp_data = rect_bivariate_interp(
        #         xin.values,
        #         yin.values,
        #         in_data[in_var].values,
        #         xout.values,
        #         yout.values,
        #     )
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
        _cmeta = {}
        # Ignore the attributes already defined for a specific variable
        # That way, individual variable metadata can override the common
        # file metadata for the dataset
        for attr in in_cfg["cmeta"]:
            if attr not in _outdata.attrs:
                _cmeta[attr] = in_cfg["cmeta"][attr]
        _outdata = _outdata.assign_attrs(_cmeta)

    return _outdata.expand_dims("time", 0)


def output_setup(
    input_file, output_file, island, out_proj, coords=None, res_m=1000.0
):
    """Create output file by copying old 1 km grid file."""
    if island == "antarctica":
        print(f"Setup Output: {input_file} -> {output_file}")
        nco = Nco()
        nco.ncks(
            input=str(input_file),
            output=str(output_file),
            options=["-x", "-v", "acab_alb,artm_alb,dzdt"],
        )
        ds_base = xr.open_dataset(output_file).load()
    else:
        # island == "greenland"
        config_file = Path("data", out_proj, f"{out_proj}grid.json")
        ds_base = build_base_greenland(
            output_file, config_file, out_proj, coords, res_m
        )

    return ds_base


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


def main(island="antarctica", resolution=1, proj_opt=None):
    """Define configuration and required variables."""
    # This mega-dictionary has all the metadata for each dataset that's going
    # to be interpolated to the output grid. Each dataset has a file,
    # some variable names, meta data for each variable, and metadata that's
    # common to all the variables in the dataset (e.g. soruce, references, etc.)
    if island == "antarctica":
        template_dset = f"{resolution}km_in"
        proj_out_name = "epsg_3031"
        projs = projections.antarctica()
        proj_out = projs[0]
    elif island == "greenland":
        template_dset = "input"
        projs = projections.greenland()
        if proj_opt is None:
            proj_out_name = "epsg_3413"
            proj_out = projs[0]
        else:
            proj_out_name = "mcb"
            proj_out = projs[1]

    config = ilib.import_module(f"{island}_config")

    input_config = config.input_config
    output_metadata = config.output_metadata
    output_variables = config.output_variables
    inout_map = config.inout_map
    ext_vars = config.ext_vars

    # Check that the files we need exist, since it takes some time to do
    # the interpolation, you hate to get halfway through and have it fail!
    all_exist = True
    for _, ds_in, in_var in inout_map:
        if not input_config[ds_in]["file"].exists():
            all_exist = False
            print(f"MISSING FILE: {input_config[ds_in]['file']}")
        else:
            print(f"FILE FOUND: {input_config[ds_in]['file']}")

    if not all_exist:
        sys.exit(1)

    # Set the output file, and coordinate variable names
    output_file = Path(
        "ncs",
        f"{island}_{resolution}km_{datetime.now().strftime('%Y_%m_%d')}.nc",
    )
    output_cfg = {"coords": {"x": "x1", "y": "y1"}}
    output = output_setup(
        input_config[template_dset]["file"],
        output_file,
        island,
        proj_out_name,
        coords=output_cfg["coords"],
    )

    # Unlink from file, since output_setup loads netCDF to memory
    # this allows writing back to the same file we open
    output.close()

    # We need both the 3412 and 3031 projections. The former is what Cryosat2
    # comes in on, the latter is what all of our output will be on
    # epsg3412, epsg3031, _ = projections.antarctica()
    grid_center_lat_lon(
        output,
        proj_out,
        proj_out_name,
        cvars=(output_cfg["coords"]["y"], output_cfg["coords"]["x"]),
    )

    if island == "antarctica":
        # Add grid information to cryosat configuration so that cryosat can
        # be transformed from its original grid to the output grid
        input_config["cryosat"]["grid_from"] = projs[1]
        input_config["cryosat"]["grid_to"] = proj_out

    # These are variables in the base file that just need metadata copied
    for dset, var in ext_vars:
        output[var] = output[var].assign_attrs(input_config[dset]["meta"][var])
        output[var] = output[var].assign_attrs(input_config[dset]["cmeta"])

    # Do the interpolation step for each variable. This could possibly be split
    # out into multiprocessing queue, since the computations are independent
    print("-" * 50)
    for out_var, ds_in, in_var in inout_map:
        output[out_var] = interp(
            input_config[ds_in], in_var, output_cfg, output
        )
        output[out_var] = output[out_var].assign_attrs(
            {"grid_mapping": proj_out_name}
        )

    # TODO: Check masking

    # The file needs metadata
    output = output.assign_attrs(output_metadata)

    # There are also coordinate variables that need metadata
    for variable in output_variables:
        if variable in ["epsg_3031", "epsg_3413", "mcb"]:
            # This creates a zero-dimension grid mapping variable so that
            # the grid mapping is cf-conventions compliant
            output[variable] = xr.Variable(
                dims=[], data=None, attrs=output_variables[variable]
            )

        # Otherwise just assign the metadata defined in the output dictionary
        if variable in output:
            output[variable] = output[variable].assign_attrs(
                output_variables[variable]
            )

    output.to_netcdf(output_file)
    print(f"-------- COMPLETE FILE SAVED TO {output_file} --------")


if __name__ == "__main__":
    main("greenland")
