#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Configuration options for Antarctic datasets.
"""
import build_antarctica as ba
from pathlib import Path
from datetime import datetime

DATA_ROOT = "/Volumes/data/piscees/ant"

input_config = {
    "1km_in": {
        "file": Path("ncs", "antarctica_1km_2017_05_03.nc"),
        "vars": ["acab_alb", "artm_alb", "dzdt"],
        "coords": {"x": "x1", "y": "y1"},
    },
    "8km_in": {
        "file": Path("ncs", "antarctica_8km_2020_03_19.nc"),
        "vars": ["acab_alb", "artm_alb", "dzdt"],
        "coords": {"x": "x1", "y": "y1"},
    },
    # Rignot Subshelf Melt rates file
    "rignot_subshelf": {
        "file": Path(
            DATA_ROOT,
            "Rignot-subShelfMeltRates/",
            "Ant_MeltingRate.flipNY.newAxes.nc",
        ),
        "load": ba.xr_load,
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
        "file": Path(
            DATA_ROOT,
            "500m.MassConsBed.AIS.Morlighem.2019",
            "BedMachineAntarctica_2019-11-05_v01.nc",
        ),
        "load": ba.xr_load,
        "vars": ["bed", "errbed", "firn", "mask", "surface", "thickness"],
        "coords": {"x": "x", "y": "y"},
        "meta": {
            "bed": {
                "long_name": "bed topography",
                "standard_name": "bedrock_altitude",
                "units": "m",
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
            "firn": {
                "long_name": "firn air content",
                "standard_name": "firn_air_content",
                "units": "m",
                "grid_mapping": "epsg_3031",
                "coordinates": "lon lat",
                "source": "REMA (Byrd Polar and Climate Research Center "
                "and the Polar Geospatial Center)",
            },
            "mask": {
                "long_name": "mask",
                "flag_values": "ocean ice_free_land grounded_ice "
                "floating_ice lake_vostok",
                "source": "Antarctic Digital Database (rock outcrop) and "
                "Jeremie Mouginot pers. comm., grounding lines)",
            },
            "surface": {
                "long_name": "ice surface elevation",
                "standard_name": "surface_altitude",
                "units": "m",
                "grid_mapping": "epsg_3031",
                "coordinates": "lon lat",
                "source": "REMA (Byrd Polar and Climate Research Center "
                "and the Polar Geospatial Center)",
            },
            "thickness": {
                "long_name": "ice thickness",
                "standard_name": "land_ice_thickness",
                "units": "m",
                "ancillary_variables": "thkerr",
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
                "Obtained from NSIDC: https://doi.org/10.5067/E1QL9HFQ7A8M "
                "Resampled from 500m grid using Nearest Neighbor -> 1km"
            ),
        },
    },
    "heatflux": {
        "file": Path(DATA_ROOT, "Martos-AIS-heatFlux", "Antarctic_GHF.xyz"),
        "load": ba.load_hf,
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
        "file": Path(
            DATA_ROOT, "Martos-AIS-heatFlux", "Antarctic_GHF_uncertainty.xyz"
        ),
        "load": ba.load_hf,
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
        "file": Path(DATA_ROOT, "cryosat2-dzdt", "CS2_dzdt.nc"),
        "load": ba.load_cryosat,
        "vars": ["dzdt", "dzdterr"],
        "coords": {"x": "x1", "y": "y1"},
        "grid_from": "epsg_3031",
        "grid_to": "epsg_3412",
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
                "comments": (
                    "Mean 1979--2010 SMB; 2D linear interpolation of 27 km "
                    "dataset."
                ),
            },
            "artm": {
                "long_name": "annual mean air temperature (2 meter)",
                "standard_name": "air_temperature",
                "units": "degree_Celsius",
                "grid_mapping": "epsg_3031",
                "coordinates": "lon lat",
                "comments": (
                    "Mean 1979--2010 t2m; 2D linear interpolation of 27 km "
                    "dataset."
                ),
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
        },
    },
    "racmo2p3": {
        "file": Path(
            DATA_ROOT,
            "RACMO2.3p2_ANT27_SMB_yearly_1979_2014",
            "RACMO2.3p2_ANT27_smb_mean_1979_2016.nc",
        ),
        "load": ba.load_racmo2p3,
        "vars": ["smb"],
        "coords": {"y": "lat", "x": "lon"},
        "meta": {
            "smb": {
                "long_name": "water equivalent surface mass balance",
                "standard_name": (
                    "land_ice_lwe_surface_specific_mass_balance"
                ),
                "units": "mm w.e. year-1",
                "grid_mapping": "epsg_3031",
                "coordinates": "lon lat",
                "comments": (
                    "Annual mean (1979-2016) gridded surface mass balance data for Antarctic"
                    "ice sheet from the regional atmospheric climate model RACMO2.3p2"
                ),
                "reference": (
                    "van Wessem, J. M., et al.: Modelling the climate and surface mass balance "
                    "of polar ice sheets using RACMO2 – Part 2: Antarctica (1979–2016), "
                    "The Cryosphere, 12, 1479–1498, https://doi.org/10.5194/tc-12-1479-2018, 2018."
                ),
            },
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
    "Mouginot-InSAR": {
        "file": Path(
            DATA_ROOT,
            "Mouginot-InSAR-450m-2019",
            "antarctic_ice_vel_phase_map_v01.nc",
        ),
        "load": ba.xr_load,
        "vars": ["VX", "VY", "ERRX", "ERRY"],
        "coords": {"x": "x", "y": "y"},
        "meta": {
            "VX": {
                "long_name": "Ice velocity in x direction",
                "standard_name": "land_ice_x_velocity",
                "units": "meter year-1",
            },
            "VY": {
                "long_name": "Ice velocity in y direction",
                "standard_name": "land_ice_y_velocity",
                "units": "meter year-1",
            },
            "ERRX": {
                "long_name": "Ice velocity in x direction error",
                "standard_name": "land_ice_x_velocity standard_error",
                "units": "meter year-1",
            },
            "ERRY": {
                "long_name": "Ice velocity in y direction error",
                "standard_name": "land_ice_y_velocity standard_error",
                "units": "meter year-1",
            },
        },
        "cmeta": {
            "source": (
                "MEaSURES Antarctica Ice Velocity Map 450m spacing\n"
                "https://nsidc.org/data/nsidc-0754/\n"
                "Department of Earth System Science, University of California, "
                "Irvine, CA, USA; Univ. Grenoble Alpes, CNRS, IRD, "
                "Grenoble INP, IGE, 38000 Grenoble, France"
            ),
            "reference": (
                "Mouginot, J., E. Rignot, and B. Scheuchl. 2019. "
                "Continent-wide, interferometric SAR phase-mapping of "
                "Antarctic ice velocity, Geophysical Research Letters. "
                "46. 9710-9718. https://doi.org/10.1029/2019GL083826"
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
            "comments": ("2D linear interpolation of provided 450 m dataset"),
        },
    },
}

# Get metadata from heatflux config (so it doesn't have to be in here twice)
input_config["heatflux_unc"]["cmeta"] = input_config["heatflux"]["cmeta"]

# Define metadata for the output file itself
output_metadata = {
    "title": "CISM-style input dataset for ice-sheet models",
    "history": "Created {} by J. Kennedy & M. Kelleher.".format(
        datetime.now().strftime("%c")
    ),
    "institution": "Oak Ridge National Laboratory",
    "references": "See https://github.com/mkstratos/cism-data",
    "Conventions": "CF-1.7",
}

# Define metadata for the output coordinate variables
output_variables = {
    "time": {
        "long_name": "time",
        "standard_name": "time",
        "axis": "T",
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
        "false_easting": 0.0,
        "false_northing": 0.0,
        "geographic_crs_name": "EPSG3031",
        "grid_mapping_name": "polar_stereographic",
        "horizontal_datum_name": "WGS84",
        "latitude_of_projection_origin": -90.0,
        "reference_ellipsoid_name": "WGS84",
        "standard_parallel": -71.0,
        "straight_vertical_longitude_from_pole": 0.0,
        "scale_factor_at_projection_origin": 1.0,
        "prime_meridian_name": "Greenwich",
        "proj4_string": (
            "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0"
            " +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
        ),
        "units": "m",
    },
}
output_cfg = {"coords": {"x": "x1", "y": "y1"}}
# Map required output variable to a dataset and input variable
# (OUTPUT VARIABLE, INPUT DATASET, INPUT VARIABLE NAME)
inout_map = [
    ("smb", "racmo2p3", "smb"),
    ("dhdt", "cryosat", "dzdt"),
    ("thk", "bedmachine", "thickness"),
    ("thkerr", "bedmachine", "errbed"),
    ("topg", "bedmachine", "bed"),
    ("bheatflx", "heatflux", "bheatflux"),
    ("bheatflxerr", "heatflux_unc", "bheatflxerr"),
    ("subm", "rignot_subshelf", "melt_actual"),
    ("subm_ss", "rignot_subshelf", "melt_steadystate"),
    ("vx", "Mouginot-InSAR", "VX"),
    ("vy", "Mouginot-InSAR", "VY"),
    ("ex", "Mouginot-InSAR", "ERRX"),
    ("ey", "Mouginot-InSAR", "ERRY"),
]

# Add metadata for extant variables (copied from the original netCDF file)
ext_vars = [
    # ("smb", "acab"),
    ("smb", "artm"),
    # ("topg", "topg"),
    # ("veloc", "verr"),
    # ("veloc", "vx"),
    # ("veloc", "vy"),
]
