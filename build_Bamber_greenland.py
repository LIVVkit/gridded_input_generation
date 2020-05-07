#!/usr/bin/env python

import os
from pathlib import Path
import datetime
import subprocess
import argparse
import xarray as xr

from util import speak
from util import finalize
from util import projections
from util.ncfunc import get_nc_file
from util import custom_argparse_types as cats

from data import bamberdem
from data import searise
from data import racmo2p3
from data import insar
from data import icebridge
from data import ice2sea
from data import csatho
from data import measures_velocity

"""
Build a CISM dataset
"""
# ==== Data Locations ====
# Link data here or edit
# ========================
DATA_ROOT = "/Volumes/data/piscees/gis"
# FIXME: Centralize data location vars so users only have to edit this info once
lc_bamber = Path(
    DATA_ROOT, "1km-res-Bamber-DEM", "Greenland_bedrock_topography_V3.nc"
)
# lc_bamber = "data/BamberDEM/Greenland_bedrock_topography_V3.nc"
lc_seaRise = "data/SeaRise/Greenland1km.nc"
lc_seaRise = Path(DATA_ROOT, "1km-res-other", "Greenland1km.nc")
# lc_racmo2p0 = "data/RACMO2.0/Racmo2MeanSMB_1961-1990.nc"
# lc_racmo2p0 = Path(
#     DATA_ROOT, "RACMO-new", "ERA-Interim", "gis_RACMO_smb_std.nc"
# )
lc_racmo2p3 = Path(
    DATA_ROOT,
    "RACMO2p3",
    "smb_rec_stats.1958-2019.BN_RACMOP2.3p2_FGRN055_GrIS.MM.nc",
)

# NOTE:  will build this file from mosaicOffsets.* files
# lc_InSAR = "data/InSAR/Joughin2015/greenland_vel_mosaic500.nc"
lc_InSAR = Path(DATA_ROOT, "joughin-InSAR-2015", "greenland_vel_mosaic500.nc")

# lc_massCon = "data/150m-MC-thickness/BedMachineGreenland-2017-09-20.nc"
lc_massCon = Path(
    DATA_ROOT, "150m-MC-thickness", "BedMachineGreenland-2017-09-20.nc"
)
# lc_massCon = "data/IceBridge/Greenland/MCdataset-2014-11-19.nc"
# lc_mask = "data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc"
lc_mask = Path(
    DATA_ROOT,
    "1km-res-other",
    "ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc",
)

# lc_csatho = (
#     "data/Csatho2014/GreenlandIceSheetdhdt_csatho/"
#     "greenland_ice_sheet_dhdt_icesat_atm_l1a_to_l2f.cdf"
# )
# lc_csatho = "data/Csatho2014/processed/Icedyndhdtave0309.nc"
lc_csatho = Path(DATA_ROOT, "dhdt-Csatho", "Icedyndhdtave0309.nc")

# ==== SETUP =====
# get args, time
# load data sets
# ================
stamp = datetime.datetime.now().strftime("%Y_%m_%d")
f_base = "templates/greenland_1km.mcb.nc"

# parse the command line arguments
parser = argparse.ArgumentParser()  # -h or --help automatically included!

parser.add_argument(
    "-s",
    "--shrink",
    type=cats.abs_existing_file,
    default="data/BamberDEM/Bambergrid_shrunk.json",
    help=(
        "JSON description of the shrunken grid specs. Use plot_grids.py "
        "to create this file."
    ),
)

volume = parser.add_mutually_exclusive_group()
volume.add_argument(
    "-v",
    "--verbose",
    help="Increase the output verbosity",
    action="store_true",
)
volume.add_argument("-q", "--quiet", help="Run silently", action="store_true")

args = parser.parse_args()

speak.notquiet(
    args, "\nBuilding the Greenland datasets in the Bamber projection."
)
speak.notquiet(
    args, "=========================================================\n"
)

# load in datasets
speak.notquiet(args, "Loading the datasets.")


nc_bamber = get_nc_file(lc_bamber, "r")
speak.verbose(args, "   Found Bamber DEM")
f_shrink = cats.abs_existing_file(args.shrink)
speak.verbose(args, "   Found shrunken Bamber grid specs")


nc_seaRise = get_nc_file(lc_seaRise, "r")
speak.verbose(args, "   Found Sea Rise data")

# nc_racmo2p0 = get_nc_file(lc_racmo2p0, "r")
# speak.verbose(args, "   Found RACMO 2.0 data")
nc_racmo2p3 = get_nc_file(lc_racmo2p3, "r")
speak.verbose(args, "   Found RACMO 2.3 data")

# try:.
#     nc_insar = get_nc_file(lc_InSAR, "r")
# except Exception:
#     speak.verbose(args, "\n   Building InSAR velocity dataset...\n")
#     subprocess.call(
#         "python util/convert_velocities.py " + os.path.dirname(lc_InSAR),
#         shell=True,
#     )
#     nc_insar = get_nc_file(lc_InSAR, "r")
# speak.verbose(args, "   Found InSAR data")


nc_massCon = get_nc_file(lc_massCon, "r")
speak.verbose(args, "   Found Mass Conserving Bed data")


nc_mask = get_nc_file(lc_mask, "r")
speak.verbose(args, "   Found Zurich mask")
nc_mask.close()

nc_csatho = get_nc_file(lc_csatho, "r")
speak.verbose(args, "   Found Csatho data")

speak.verbose(args, "\n   All data files found!")


# ===== Bamber DEM ======
# this is a 1km dataset
# =======================
speak.verbose(args, "\nBuilding the base dataset: " + f_base)

speak.notquiet(args, "\nCreating the base grid."),

nc_base, base = bamberdem.build_base(f_base, nc_bamber)

speak.notquiet(args, "   Done!")

# ==== Projections ====
# All the projections
# needed for the data
# =====================
speak.notquiet(args, "\nGetting the projections.")

proj_epsg3413, proj_eigen_gl04c = projections.greenland()

speak.notquiet(args, "   Done!")

# transform meshes.
speak.verbose(
    args, "   Creating the transform meshes: base Bamber grid to EPSG-3413."
)

trans = projections.transform(base, proj_eigen_gl04c, proj_epsg3413)

speak.notquiet(args, "   Done!")

speak.notquiet(args, "\nAdding the lon / lat dimension")
projections.grid_center_latlons(
    nc_base, base, proj_eigen_gl04c, "mcb", ("y", "x")
)
speak.notquiet(args, "   Done!")

# ==== SeaRise Data =====
# this is a 1km dataset
# =======================
speak.notquiet(args, "\nGetting bheatflx and presartm from the SeaRise data.")

searise.bheatflx_artm_bamber(args, nc_seaRise, nc_base, base)

nc_seaRise.close()
# ==== RACMO2.0 Data =====
# this is a 1km dataset
# ========================
speak.notquiet(args, "\nGetting acab from the RACMO 2.3 data.")
racmo2p3.acab_bamber(args, nc_racmo2p3, nc_base, base)
nc_racmo2p3.close()
# ==== InSAR velocity Data ====
# this is a 500m dataset in
# the ESPG-3413 projection
# =============================
speak.notquiet(args, "\nGetting vy, vx, ey, and ex from the InSAR data.")

measures_velocity.velocity(args, nc_base, base, proj_eigen_gl04c)

# nc_insar.close()
# ==== Mass Conserving Bed Data ===
# This is the new (2015) bed data
# =================================
speak.notquiet(
    args,
    "\nGetting thk, topg, and topgerr from the mass conserving bed data.",
)

icebridge.mcb_bamber(
    args,
    nc_massCon,
    nc_bamber,
    nc_base,
    base,
    trans,
    proj_eigen_gl04c,
    proj_epsg3413,
)

nc_bamber.close()
nc_massCon.close()

# == Csatho dHdt data needs to follow thickness interpolation) ==
# Approx. 8km dataset
# =====================
speak.notquiet(args, "\nGetting dhdt from the Csatho data.")
csatho.dhdt_all(
    args,
    nc_csatho,
    nc_base,
    base,
    proj_out=None,  # proj_epsg3413, proj_eigen_gl04c
)
speak.notquiet(args, "   Done!")
nc_csatho.close()

# ==== Zurich mask =====
# apply mask, and get
# new surface variable
# ======================
speak.notquiet(args, "\nGetting the Zurich Mask.")

base = None
nc_base.close()  # need to read in some data from nc_base now
nc_base = get_nc_file(f_base, "r+")

nc_mask = xr.open_dataset(lc_mask)
ice2sea.apply_mask(args, nc_mask, nc_base)

nc_mask.close()
# ==== Done getting data ====
# ===========================
for var in nc_base.variables:
    if var not in ["mcb", "lon", "lat", "x1", "y1"]:
        nc_base[var].grid_mapping = "mcb"

nc_base.close()

speak.notquiet(args, "   Done!")

# ==== add time dim and shrink ====
# apply to all the variables and
# shrink to size around ice sheet
# =================================
speak.notquiet(
    args, "\nAdding the time dimension and creating the 1km dataset."
)

f_1km = "complete/greenland_1km_" + stamp + ".mcb.nc"
f_template = "templates/greenland.mcb.config"
# Write file-level metadata
output_metadata = {
    "title": "CISM-style input dataset for ice-sheet models",
    "history": "Created {} by M. Kelleher & J. Kennedy.".format(
        datetime.datetime.now().strftime("%Y-%m-%d")
    ),
    "institution": "Oak Ridge National Laboratory",
    "references": "https://github.com/LIVVkit/gridded_input_generation",
    "Conventions": "CF-1.7",
}

finalize.add_time_and_shrink(
    args, "mcb", f_base, f_1km, f_template, f_shrink, output_metadata
)

# ==== Coarsen ====
# make 2, 4 and 8
# km datasets
# ==================
speak.notquiet(args, "\nCreating coarser datasets.")

coarse_list = [2, 4, 5, 8]  # in km

finalize.coarsen(args, "mcb", f_1km, f_template, coarse_list)

# ==== and done! ====
# ===================
speak.notquiet(args, "\nFinished building the datasets.")
