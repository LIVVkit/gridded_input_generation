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
    "seaRise": {
        "file": Path(DATA_ROOT, "SeaRise", "Greenland1km.nc"),
        "load": ba.load_cryosat,
        "vars": ["bheatflx", "presartm", "presprcp", "thk", "topg"],
        "coords": {"x": "x", "y": "y"},
        "grid_from": "mcb",
        "grid_to": "mcb",
        "meta": {
            "bheatflx": {
                "long_name": "Basal Heat Flux",
                "standard_name": (
                    "upward_geothermal_heat_flux_at_ground_level_in_land_ice"
                ),
                "units": "W m-2",
                "reference": (
                    "Shapiro, N.M. and Ritzwoller, M.H. (2004), Inferring "
                    "Surface "
                    "Heat Flux Distributions Guided by a Global Seismic Model: "
                    "Particular Application to Antarctica; Earth and Planetary "
                    "Science letters, 223: 213-224."
                ),
                "source": (
                    "http://ciei.colorado.edu/~nshapiro/MODEL/ASC_VERSION/index.html"
                ),
            },
            "presartm": {
                "reference": (
                    "Ettema J., M.R. van den Broeke, E. van Meigaard, "
                    "W.J. van de Berg, J.L. Bamber, J.E. Box, and R.C. Bales "
                    "(2009), Higher surface mass balance of the Greenland ice "
                    "sheet revealed by high-resolution climate modeling, "
                    "Geophys. Res. Lett., 36, L12501, doi:10.1029/2009GL038110"
                ),
                "source": (
                    "Janneke Ettema; Institute for Marine and Atmospheric "
                    "Research, Utrecht University, Utrecht, Netherlands"
                ),
                "long_name": "Annual Mean Air Temperature (2 meter)",
                "standard_name": "air_temperature",
                "units": "degree_Celsius",
            },
            "presprcp": {
                "reference": (
                    "Ettema J., M.R. van den Broeke, E. van Meigaard, "
                    "W.J. van de Berg, J.L. Bamber, J.E. Box, and R.C. Bales "
                    "(2009), Higher surface mass balance of the Greenland ice "
                    "sheet revealed by high-resolution climate modeling, "
                    "Geophys. Res. Lett., 36, L12501, doi:10.1029/2009GL038110"
                ),
                "source": (
                    "Janneke Ettema; Institute for Marine and Atmospheric "
                    "Research, Utrecht University, Utrecht, Netherlands"
                ),
                "long_name": "Present Precipitation",
                "standard_name": "lwe_precipitation_rate",
                "units": "m year-1",
            },
            "thk": {
                "reference": (
                    "[1a] Bamber, J.L., R.L. Layberry, S.P. Gogenini. 2001. "
                    "A new ice thickness and bed data set for the Greenland "
                    "ice sheet 1: Measurement, data reduction, and errors. "
                    "Journal of Geophysical Research 106 (D24): 33773-33780. "
                    "[1b] Bamber, J.L., R.L. Layberry, S.P. Gogenini. 2001. "
                    "A new ice thickness and bed data set for the Greenland "
                    "ice sheet 2: Relationship between dynamics and basal "
                    "topography. Journal of Geophysical Research 106 (D24): "
                    "33781-33788."
                ),
                "source": (
                    "[1] http://nsidc.org/data/nsidc-0092.html, "
                    "[2] https://www.cresis.ku.edu/data/Greenland"
                ),
                "comments": (
                    "The ice thickness field incorporates (corrected) 5 km "
                    "gridded data from the National Snow and Ice Data Center "
                    "(NSIDC), gridded data for the Jakobshavn and Pettermann "
                    "Glacier regions from the Center for Remote Sensing of "
                    "Ice Sheets (CReSIS), and flightline data for the "
                    "Kangerdlugsuaq and Helheim Glacier regions also "
                    "from CReSIS.  Additional data points with zero thickness "
                    "were introduced to smooth the ice margin but artifacts "
                    "of the 5 km grid remain in some locations.  Because of "
                    "the large number of data points, interpolation onto the "
                    "output grid (using natgridd's nonlinear natural neighbor "
                    "algorithm) was performed in blocks to reduce the "
                    "computing time required.  The available data for the "
                    "Jakobshavn Glacier region was for bed elevation only; "
                    "To calculate thicknesses from these bed elevations, "
                    "the bed elevations were subtracted from surface "
                    "elevations obtained using a bicubic spline of 5km gridded "
                    "surface data from NSIDC."
                ),
                "long_name": "Ice Thickness",
                "standard_name": "land_ice_thickness",
                "units": "meters",
            },
            "topg": {
                "reference": (
                    "[1a] Bamber, J.L., R.L. Layberry, S.P. Gogenini. 2001. "
                    "A new ice thickness and bed data set for the Greenland "
                    "ice sheet 1: Measurement, data reduction, and errors. "
                    "Journal of Geophysical Research 106 (D24): 33773-33780. "
                    "[1b] Bamber, J.L., R.L. Layberry, S.P. Gogenini. 2001. "
                    "A new ice thickness and bed data set for the Greenland "
                    "ice sheet 2: Relationship between dynamics and basal "
                    "topography. Journal of Geophysical Research 106 (D24): "
                    "33781-33788."
                ),
                "source": (
                    "[1] http://nsidc.org/data/nsidc-0092.html, "
                    "[2] http://www.gebco.net/data_and_products/gridded_bathymetry_data/, "
                    "[3] https://www.cresis.ku.edu/data/Greenland"
                ),
                "comments": (
                    "The bedrock topography field incorporates (corrected) "
                    "5 km gridded data from the National Snow and Ice Data "
                    "Center (NSIDC), bathymetric data and surface elevation "
                    "data for Ellesmere Island from the General Bathymetric "
                    "Chart of the Oceans (GEBCO), gridded data for the "
                    "Jakobshavn and Pettermann Glacier regions from the Center "
                    "for Remote Sensing of Ice Sheets (CReSIS), and flightline "
                    "data for the Kangerdlugsuaq and Helheim Glacier regions "
                    "also from CReSIS.  Because of the large number of data "
                    "points, interpolation onto the output grid (using "
                    "natgridd's nonlinear natural neighbor algorithm) was "
                    "performed in blocks to reduce the computing time "
                    "required.  The available flightline data was for ice "
                    "thickness only; To calculate bed elevations from these "
                    "thicknesses, the thicknesses were subtracted from "
                    "surface elevations obtained using a bicubic spline of "
                    "5km gridded surface data from NSIDC."
                ),
                "long_name": "Bedrock Topography",
                "standard_name": "bedrock_altitude",
                "units": "meters",
            },
        },
        "cmeta": {},
    },
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
        "load": ba.load_cryosat,
        "grid_from": "epsg_3413",
        "grid_to": "mcb",
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

output_metadata = {
    "title": "CISM-style input dataset for ice-sheet models",
    "history": "Created {} by J. Kennedy & M. Kelleher.".format(
        datetime.now().strftime("%c")
    ),
    "institution": "Oak Ridge National Laboratory",
    "references": "See https://github.com/mkstratos/cism-data",
    "Conventions": "CF-1.7",
}

output_variables = {
    "time": {
        "long_name": "time",
        "standard_name": "time",
        "axis": "T",
        "units": "common_years since 2009-01-01 00:00:00",
        "calendar": "365_day",
        "comments": (
            "The initial time here is an estimate of the nominal date for "
            "Joughin's 2015 InSAR data. Because this is a synthesis of "
            "datasets across many time periods, the inital date is inherently "
            "fuzzy and should be changed to suit your purposes."
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
    "epsg_3413": {
        "false_easting": 0.0,
        "false_northing": 0.0,
        "geographic_crs_name": "EPSG3413",
        "grid_mapping_name": "polar_stereographic",
        "horizontal_datum_name": "WGS84",
        "latitude_of_projection_origin": 90.0,
        "reference_ellipsoid_name": "WGS84",
        "standard_parallel": 70.0,
        "straight_vertical_longitude_from_pole": -45.0,
        "scale_factor_at_projection_origin": 1.0,
        "prime_meridian_name": "Greenwich",
        "proj4_string": (
            "+proj=stere +lat_ts=70.0 +lat_0=90 +lon_0=-45.0 +k_0=1.0 +x_0=0.0 "
            "+y_0=0.0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
        ),
        "units": "m",
    },
    "mcb": {
        "false_easting": 0.0,
        "false_northing": 0.0,
        "geographic_crs_name": "EPSG3413",
        "grid_mapping_name": "polar_stereographic",
        "horizontal_datum_name": "WGS84",
        "latitude_of_projection_origin": 90.0,
        "reference_ellipsoid_name": "WGS84",
        "standard_parallel": 71.0,
        "straight_vertical_longitude_from_pole": 321.0,
        "scale_factor_at_projection_origin": 1.0,
        "prime_meridian_name": "Greenwich",
        "proj4_string": (
            "+proj=stere +lat_ts=71.0 +lat_0=90 +lon_0=321.0 +k_0=1.0"
        ),
        "units": "m",
    },
}
output_cfg = {"coords": {"x": "x1", "y": "y1"}}
inout_map = [
    ("thk", "massCon", "thickness"),
    ("bheatflx", "seaRise", "bheatflx"),
    ("artm", "seaRise", "presartm"),
]
ext_vars = []
