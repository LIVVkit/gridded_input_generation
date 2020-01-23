#!/usr/bin/env python2

"""
Build a CISM dataset
"""

import os
import datetime
import subprocess
import argparse

from util import speak
from util import finalize
from util import projections
from util.ncfunc import get_nc_file
from util import custom_argparse_types as cats

# from data import bamberdem
# from data import ice2sea
from data import epsg3413
from data import searise
from data import racmo2p3
from data import insar
from data import icebridge
from data import csatho

# === Data Locations ====
# Link data here or edit 
# =======================
# FIXME: Centralize data location vars so users only have to edit this info once
lc_epsg = 'data/EPSG3413/EPSG3413grid.json'  # NOTE: created by plot_grids.py
lc_bamber = 'data/BamberDEM/Greenland_bedrock_topography_V3.nc'
lc_seaRise = 'data/SeaRise/Greenland1km.nc'
lc_racmo2p3 = 'data/RACMO2.3/smb_1km_GrIS_downscaled_RACMO2_3.nc'
lc_InSAR = 'data/InSAR/Joughin2015/greenland_vel_mosaic500.nc'  # NOTE:  will build this file from mosaicOffsets.* files
lc_massCon = 'data/IceBridge/Greenland/MCdataset-2014-11-19.nc'
lc_mask = 'data/Ice2Sea/ice2sea_Greenland_geometry_icesheet_mask_Zurich.nc'
lc_csatho = 'data/Csatho2014/GreenlandIceSheetdhdt_csatho/greenland_ice_sheet_dhdt_icesat_atm_l1a_to_l2f.cdf'

# === SETUP ====
# get args, time
# load data sets
# ==============
stamp = datetime.date.today().strftime("%Y_%m_%d")
f_1km = 'greenland_1km_'+stamp+'.epsg3413.nc'
f_base = 'templates/greenland_1km.epsg3413.nc'
f_template = 'templates/greenland.epsg3413.config'

# parse the command line arguments
parser = argparse.ArgumentParser()  # -h or --help automatically included!

parser.add_argument('-s', '--shrink', 
                    type=cats.abs_existing_file,
                    default='data/EPSG3413/EPSG3413grid_shrunk.json',
                    help='JSON description of the shrunken grid specs. Use '
                         'plot_grids.py to create this file.')

parser.add_argument('-o', '--out-dir',
                    type=cats.abs_creation_dir,
                    default=os.path.join(os.getcwd(), 'complete'),
                    help='Location to output complete datasets.')

parser.add_argument('-t', '--use-template', action='store_true', 
                    help='Reuse the current data template.')

volume = parser.add_mutually_exclusive_group()
volume.add_argument('-v', '--verbose', action='store_true',
                    help='Increase the output verbosity')
volume.add_argument('-q', '--quiet', action='store_true',
                    help='Run silently')

args = parser.parse_args()

f_1km = os.path.join(args.out_dir, f_1km)

speak.notquiet(args, '\nBuilding the Greenland datasets in the EPSG:3413 projection.')
speak.notquiet(args,  '============================================================\n')

# load in datasets
speak.notquiet(args, 'Loading the datasets.')

f_epsg = cats.abs_existing_file(lc_epsg)
speak.verbose(args, '   Found EPSG:3413 grid specs')
f_shrink = cats.abs_existing_file(args.shrink)
speak.verbose(args, '   Found shrunken EPSG:3413 grid specs')

nc_bamber = get_nc_file(lc_bamber, 'r')
speak.verbose(args, '   Found Bamber DEM')
    
if not args.use_template:
    nc_seaRise = get_nc_file(lc_seaRise, 'r')
    speak.verbose(args, '   Found Sea Rise data')

    nc_racmo2p3 = get_nc_file(lc_racmo2p3, 'r')
    speak.verbose(args, '   Found RACMO 2.3 data')

    try:
        nc_insar = get_nc_file(lc_InSAR, 'r')
    except IOError:
        speak.verbose(args, '\n   Building InSAR velocity dataset...\n')
        subprocess.call('python util/convert_velocities.py '+os.path.dirname(lc_InSAR), shell=True)
        nc_insar = get_nc_file(lc_InSAR, 'r')
    speak.verbose(args, '   Found InSAR data')

    nc_massCon = get_nc_file(lc_massCon, 'r')
    speak.verbose(args, '   Found Mass Conserving Bed data')

    nc_mask = get_nc_file(lc_mask, 'r')
    speak.verbose(args, '   Found Zurich mask')

    nc_csatho = get_nc_file(lc_csatho, 'r')
    speak.verbose(args, '   Found Csatho data')

    speak.verbose(args, '\n   All data files found!')

    # ===== EPSG:3413 ======
    # this is a 1km dataset 
    # ======================
    speak.notquiet(args, '\nCreating the base grid.'),
    nc_base, base = epsg3413.build_base(f_base, f_epsg, 1000.)
    speak.notquiet(args, '   Done!')

    # === Projections ====
    # All the projections 
    # needed for the data 
    # ====================
    speak.notquiet(args, '\nGetting the projections.')
    proj_epsg3413, proj_eigen_gl04c = projections.greenland()
    speak.notquiet(args, '   Done!')

    # ===== Lat,Lon ======
    # Get the cell-center 
    # lats and lons       
    # ====================
    speak.notquiet(args, '\nDetermining the  latitudes and longitudes of the grid cell-centers.')
    projections.grid_center_latlons(nc_base, base, proj_epsg3413, 'epsg_3413')
    speak.notquiet(args, '   Done!')

    # === SeaRise Data =====
    # this is a 1km dataset
    # ======================
    speak.notquiet(args, '\nGetting bheatflx and presartm from the SeaRise data.')
    searise.bheatflx_artm_epsg3413(args, nc_seaRise, nc_base, base, proj_epsg3413, proj_eigen_gl04c)
    speak.notquiet(args, '   Done!')
    nc_seaRise.close()

    # === RACMO2.3 Data =====
    # this is a 1km dataset
    # =======================
    speak.notquiet(args, '\nGetting acab from the RACMO 2.0 data.')
    racmo2p3.acab_epsg3413(args, nc_racmo2p3, nc_base, base)
    speak.notquiet(args, '   Done!')
    nc_racmo2p3.close()

    # === InSAR velocity Data ====
    # this is a 500m dataset in
    # the ESPG-3413 projection
    # ============================
    speak.notquiet(args, '\nGetting vy, vx, ey, and ex from the InSAR data.')
    insar.velocity_epsg3413(args, nc_insar, nc_base, base)
    speak.notquiet(args, '   Done!')
    nc_insar.close()

    # === Mass Conserving Bed Data ===
    # This is the new (2015) bed data
    # ================================
    speak.notquiet(args, '\nGetting thk, topg, and topgerr from the mass conserving bed data.')
    icebridge.mcb_epsg3413(args, nc_massCon, nc_bamber, nc_base, base, proj_epsg3413, proj_eigen_gl04c)
    speak.notquiet(args, '   Done!')
    nc_bamber.close()
    nc_massCon.close()

    # = Csatho dHdt data ==
    # Approx. 8km dataset
    # =====================
    speak.notquiet(args, '\nGetting dhdt from the Csatho data.')
    csatho.dhdt_epsg3413(args, nc_csatho, nc_base, base, proj_epsg3413)
    speak.notquiet(args, '   Done!')
    nc_csatho.close()

    # === Done getting data ====
    # ==========================
    nc_base.close()

    # === Generate Mask ====
    # ======================
    speak.notquiet(args, '\nBuilding mask.')
    import segment_topo
    # NOTE: the space in front of altitude so argparse doesn't think it's an option
    seg_arg_list = ['-w', '-i', f_base, '-a', ' -200.']
    if args.quiet:
        seg_arg_list.append('-q')
    elif args.verbose:
        seg_arg_list.append('-v')
    seg_args = segment_topo.parse_args(seg_arg_list)
    segment_topo.main(seg_args)
    speak.notquiet(args, '   Done!')

# === add time dim and shrink ====
# apply to all the variables and  
# shrink to size around ice sheet 
# ================================
speak.notquiet(args, '\nAdding the time dimension and creating the 1km dataset.')
finalize.add_time_and_shrink(args, 'epsg_3413', f_base, f_1km, f_template, f_shrink)

# === Coarsen ==== 
# make 2, 4 and 8  
# km datasets      
# =================
speak.notquiet(args, '\nCreating coarser datasets.')
coarse_list = [2, 4, 5, 8]   # in km
finalize.coarsen(args, 'epsg_3413', f_1km, f_template, coarse_list)

# === and done! ====
# ==================
speak.notquiet(args, '\nFinished building the datasets.')
