#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Configuration options for Greenland datasets.
"""
import build_antarctica as ba
from pathlib import Path
from datetime import datetime

DATA_ROOT = "data"

input_config = {
    "input": {"file": None},
    "bamber": {
        "file": Path(
            DATA_ROOT, "BamberDEM", "Greenland_bedrock_topography_V3.nc"
        ),
    },
    "seaRise": {"file": Path(DATA_ROOT, "SeaRise", "Greenland1km.nc")},
    "racmo2p0": {
        "file": Path(DATA_ROOT, "RACMO2.3", "Racmo2MeanSMB_1961-1990.nc")
    },
    # NOTE:  will build this file from mosaicOffsets.* files
    "InSAR": {
        "file": Path(
            DATA_ROOT, "InSAR", "Joughin2015", "greenland_vel_mosaic500.nc"
        )
    },
    "massCon": {
        "file": Path(
            DATA_ROOT, "150m-MC-thickness", "BedMachineGreenland-2017-09-20.nc"
        ),
        "load": ba.xr_load,
        "vars": ["thickness", "bed", "mask", "errbed"],
        "coords": {"x": "x", "y": "y"},
        "meta": {
            "thickness": {
                "long_name": "ice thickness",
                "standard_name": "land_ice_thickness",
                "units": "meters",
                "grid_mapping": "mapping",
                "source": "Mass conservation (Mathieu Morlighem)",
            }
        },
    },
    "mask": {
        "file": Path(
            DATA_ROOT,
            "Ice2Sea" "ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc",
        )
    },
}
output_metadata = {}
output_variables = {}
inout_map = [("thk", "massCon", "thickness")]
ext_vars = []
