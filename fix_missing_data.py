#!/usr/bin/env python

import glob
from datetime import datetime
import xarray as xr
import numpy as np

date = "2020_02_13"
grid = "epsg3413"

data_files = sorted(glob.glob(f"complete/greenland_*km_{date}.{grid}.nc"))

masked_vars_all = {
    "mcb": {
        "acab": 1e10,
        "ex": 1e10,
        "ey": 1e10,
        "topg": -99,
        "topgerr": 100.0,
        "usrf": -99,
        "vx": 0.0,
        "vy": 0.0,
    },
    "epsg3413": {
        "artm": 0.0,
        "ex": 1e10,
        "ey": 1e10,
        "topg": -99,
        "vx": 0.0,
        "vy": 0.0,
        "dhdt": 0,
    },
}
masked_vars = masked_vars_all[grid]

for data_file in data_files:
    dat = xr.open_dataset(data_file, decode_times=False).load()
    dat.close()
    print(dat)
    for idx, var in enumerate(masked_vars):
        dat[var] = dat[var].where(dat[var] < 1e30)
        if var == "topgerr":
            dat[var] = dat[var].where(dat[var] != 0)
        dat[var] = dat[var].fillna(masked_vars[var])
        dat[var] = dat[var].assign_attrs(
            {
                "missing_value": np.float32(masked_vars[var]),
                "_FillValue": np.float32(masked_vars[var]),
            }
        )
dat["x1"] = dat["x1"].assign_attrs({"axis": "X"})
dat["y1"] = dat["y1"].assign_attrs({"axis": "Y"})

if grid == "mcb":
    dat["usrf"] = dat["usrf"].assign_attrs({"units": "m"})
    dat["acab"] = dat["acab"].assign_attrs(
        {"_FillValue": dat["acab"].missing_value}
    )
    dat["bheatflx"] = dat["bheatflx"].assign_attrs(
        {
            "standard_name": "upward_geothermal_heat_flux_at_ground_level_in_land_ice"
        }
    )

if grid == "epsg3413":
    _inmap = "epsg_3413"
    mapping_var = "epsg_3413"
else:
    _inmap = "mcb"
    mapping_var = "mapping"

attrs = dat[_inmap].attrs
dat = dat.drop(_inmap)
dat[mapping_var] = xr.Variable(data=None, dims=[], attrs=attrs)

for var in dat.data_vars:
    dat[var] = dat[var].assign_attrs({"grid_mapping": mapping_var})

for var in ["x1", "y1", "lon", "lat"]:
    try:
        del dat[var].attrs["grid_mapping"]
    except KeyError:
        continue

output_metadata = {
    "title": "CISM-style input dataset for ice-sheet models",
    "history": "Created {} by J. Kennedy & M. Kelleher.".format(
        datetime.now().strftime("%c")
    ),
    "institution": "Oak Ridge National Laboratory",
    "references": "See https://github.com/mkstratos/cism-data",
    "Conventions": "CF-1.7",
}
dat = dat.assign_attrs(output_metadata)
dat.to_netcdf(data_file)
