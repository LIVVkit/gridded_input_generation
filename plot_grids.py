#!/usr/bin/env python

"""
Plot the Bamber and EPSG grids. 
"""

import json
import argparse

import numpy as np
from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot as plt

from util import speak
from util import projections
from util.ncfunc import get_nc_file

# ==== Data Locations ====
# Link data here or edit
# ========================
# FIXME: Centralize data location vars so users only have to edit this info once
lc_template = "data/BamberDEM/Greenland_bedrock_topography_V3.nc"
# lc_base = "complete/mcb_datasets/greenland_1km_2016_12_01.mcb.nc"
lc_base = "complete/greenland_1km_2020_02_13.epsg3413.nc"
# ========================


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-s", "--show", help="Show the generated plot.", action="store_true"
    )

    volume = parser.add_mutually_exclusive_group()
    volume.add_argument(
        "-v",
        "--verbose",
        help="Increase the output verbosity",
        action="store_true",
    )
    volume.add_argument(
        "-q", "--quiet", help="Run silently", action="store_true"
    )

    return parser.parse_args(args)


def write_grid_spec(file_name, ll, ur, xs, ys, lons, lats, proj):
    grid = {
        "ll": ll,
        "ur": ur,
        "xs": xs,
        "ys": ys,
        "lons": lons,
        "lats": lats,
        "projstring": projections.proj_string(proj),
    }

    for key, val in grid.items():
        if isinstance(val, np.ndarray):
            grid[key] = [float(ival) for ival in val]
        elif not isinstance(val, list) and key != "projstring":
            grid[key] = [float(ival) for ival in val]

    with open(file_name, "w") as f:
        json.dump(grid, f, indent=2)


def main(args):
    speak.notquiet(
        args,
        "\nPlotting the representation of the Bamber grid in the EPSG:3413 projection",
    )
    speak.notquiet(
        args,
        "==========================================================================\n",
    )

    speak.notquiet(args, "Loading the datasets.")
    nc_template = get_nc_file(lc_template, "r")
    nc_base = get_nc_file(lc_base, "r")
    speak.verbose(args, "   Found Bamber DEM")
    speak.verbose(args, "\n   All data files found!")

    speak.notquiet(args, "\nCreating the 1 km Bamber template grid."),
    template = projections.DataGrid()
    template.y = nc_template.variables["projection_y_coordinate"]
    template.x = nc_template.variables["projection_x_coordinate"]
    template.ny = template.y[:].shape[0]
    template.nx = template.x[:].shape[0]
    template.make_grid()
    speak.notquiet(args, "   Done!")

    speak.notquiet(args, "\nCreating the 1 km Bamber grid."),
    base = projections.DataGrid()
    base.y = nc_base.variables["y1"]
    base.x = nc_base.variables["x1"]
    base.ny = base.y[:].shape[0]
    base.nx = base.x[:].shape[0]
    base.make_grid()
    speak.notquiet(args, "   Done!")

    speak.notquiet(args, "\nGetting the projections.")
    proj_epsg3413, proj_eigen_gl04c = projections.greenland()
    speak.notquiet(args, "   Done!")

    speak.notquiet(args, "\nProject the Bamber template grid into EPSG:3413.")

    T_ll = (np.amin(template.x_grid[0, :]), np.amin(template.y_grid[0, :]))
    T_lr = (np.amax(template.x_grid[0, :]), np.amax(template.y_grid[0, :]))
    T_ur = (np.amax(template.x_grid[-1, :]), np.amax(template.y_grid[-1, :]))
    T_ul = (np.amin(template.x_grid[-1, :]), np.amin(template.y_grid[-1, :]))

    T_xs = np.array([T_ll[0], T_lr[0], T_ur[0], T_ul[0], T_ll[0]])
    T_ys = np.array([T_ll[1], T_lr[1], T_ur[1], T_ul[1], T_ll[1]])

    T_lons, T_lats = proj_eigen_gl04c(T_xs, T_ys, inverse=True)

    T2E_xs, T2E_ys = proj_epsg3413(T_lons, T_lats)

    e_ll = (
        np.floor(np.mean([T2E_xs[0], T2E_xs[3]]) / 1000.0) * 1000.0,
        np.floor(np.mean([T2E_ys[0], T2E_ys[1]]) / 1000.0) * 1000.0,
    )

    e_ur = (
        np.ceil(np.mean([T2E_xs[1], T2E_xs[2]]) / 1000.0) * 1000.0,
        np.ceil(np.mean([T2E_ys[3], T2E_ys[2]]) / 1000.0) * 1000.0,
    )

    TE_xs = np.array([e_ll[0], e_ur[0], e_ur[0], e_ll[0], e_ll[0]])
    TE_ys = np.array([e_ll[1], e_ll[1], e_ur[1], e_ur[1], e_ll[1]])

    TE_lons, TE_lats = proj_epsg3413(TE_xs, TE_ys, inverse=True)

    write_grid_spec(
        "Bambergrid.json",
        T_ll,
        T_ur,
        T_xs,
        T_ys,
        T_lons,
        T_lats,
        proj_eigen_gl04c,
    )
    write_grid_spec(
        "EPSG3413grid.json",
        [TE_xs[0], TE_ys[0]],
        [TE_xs[2], TE_ys[2]],
        TE_xs,
        TE_ys,
        TE_lons,
        TE_lats,
        proj_epsg3413,
    )

    speak.notquiet(args, "\nProject the Bamber grid into EPSG:3413.")

    B_ll = (np.amin(base.x_grid[0, :]), np.amin(base.y_grid[0, :]))
    B_lr = (np.amax(base.x_grid[0, :]), np.amax(base.y_grid[0, :]))
    B_ur = (np.amax(base.x_grid[-1, :]), np.amax(base.y_grid[-1, :]))
    B_ul = (np.amin(base.x_grid[-1, :]), np.amin(base.y_grid[-1, :]))

    B_xs = np.array([B_ll[0], B_lr[0], B_ur[0], B_ul[0], B_ll[0]])
    B_ys = np.array([B_ll[1], B_lr[1], B_ur[1], B_ul[1], B_ll[1]])

    B_lons, B_lats = proj_eigen_gl04c(B_xs, B_ys, inverse=True)

    B2E_xs, B2E_ys = proj_epsg3413(B_lons, B_lats)
    _locs = (
        ("Lower Left", 0),
        ("Lower Right", 1),
        ("Upper Right", 2),
        ("Upper Left", 3),
    )
    speak.verbose(args, "    Bamber in EPSG:3413:")
    for desc, idx in _locs:
        speak.verbose(
            args, f"    {desc} (x,y): ({B2E_xs[idx]}, {B2E_ys[idx]})\n"
        )

    speak.notquiet(
        args,
        "\nDeterniming the new EPSG:3413 grid from the transformed Bamber grid.",
    )

    e_ll = [
        np.floor(np.mean([B2E_xs[0], B2E_xs[3]]) / 1000.0) * 1000.0,
        np.floor(np.mean([B2E_ys[0], B2E_ys[1]]) / 1000.0) * 1000.0,
    ]

    e_ur = [
        np.ceil(np.mean([B2E_xs[1], B2E_xs[2]]) / 1000.0) * 1000.0,
        np.ceil(np.mean([B2E_ys[3], B2E_ys[2]]) / 1000.0) * 1000.0,
    ]

    E_xs = np.array([e_ll[0], e_ur[0], e_ur[0], e_ll[0], e_ll[0]])
    E_ys = np.array([e_ll[1], e_ll[1], e_ur[1], e_ur[1], e_ll[1]])

    E_lons, E_lats = proj_epsg3413(E_xs, E_ys, inverse=True)

    # NOTE: Best normal grid
    offset_y128 = -19000.0
    # NOTE: Working extended grid
    # offset_y128 = -101000.0
    # NOTE: Idea extended grid
    # offset_y128 = -84000.0

    # NOTE: Best normal grid
    offset_xl = -80000.0
    offset_x128 = -48000.0
    # NOTE: Best extended grid
    # offset_xl = -366000.0
    # offset_x128 =-130000.0
    # NOTE: Ideal extended grid
    # offset_xl = -400000.0
    # offset_x128 =-130000.0

    O_ys = np.array(
        [
            e_ll[1] + offset_y128,
            e_ll[1] + offset_y128,
            e_ur[1] - offset_y128,
            e_ur[1] - offset_y128,
            e_ll[1] + offset_y128,
        ]
    )
    O_xs = np.array(
        [
            e_ll[0] + offset_xl + offset_x128,
            e_ur[0] - offset_x128 - 1000.0,
            e_ur[0] - offset_x128 - 1000.0,
            e_ll[0] + offset_xl + offset_x128,
            e_ll[0] + offset_xl + offset_x128,
        ]
    )

    O_lons, O_lats = proj_epsg3413(O_xs, O_ys, inverse=True)
    speak.notquiet(args, "    New EPSG:3413 grid:")

    for desc, idx in _locs:
        if idx != 0:
            _strout = "\n"
        else:
            _strout = ""
        _strout += f"    {desc} (x, y): ({O_xs[idx]}, {O_ys[0]})\n"
        _strout += f"    {desc} (lat, lon): ({O_lats[idx]}, {O_lons[0]})"
        speak.notquiet(args, _strout)

    grid_points_y = int((O_ys[3] - O_ys[0]) / 1000.0 + 1)
    grid_points_x = int((O_xs[1] - O_xs[0]) / 1000.0 + 1)
    speak.notquiet(
        args, "\n      Number of (1km) grid points in y: " + str(grid_points_y)
    )
    speak.notquiet(
        args, "      Number of (1km) grid points in x: " + str(grid_points_x)
    )

    grid_p128_y = grid_points_y + (grid_points_y % 128)
    grid_p128_x = grid_points_x + (grid_points_x % 128)
    speak.notquiet(
        args,
        "\n      Number of (1km) grid points for next multiple of 128 in y: "
        + str(grid_p128_y),
    )
    speak.notquiet(
        args,
        "      Number of (1km) grid points for next multiple of 128 in x: "
        + str(grid_p128_x),
    )
    speak.notquiet(
        args,
        "      Number of km needed to (add, subtract) to y: "
        + str(128 - (grid_points_y % 128))
        + ", "
        + str(grid_points_y % 128),
    )
    speak.notquiet(
        args,
        "      Number of km needed to (add, subtract) to x: "
        + str(128 - (grid_points_x % 128))
        + ", "
        + str(grid_points_x % 128),
    )

    write_grid_spec(
        "Bambergrid_shrunk.json",
        B_ll,
        B_ur,
        B_xs,
        B_ys,
        B_lons,
        B_lats,
        proj_eigen_gl04c,
    )
    write_grid_spec(
        "EPSG3413grid_shrunk.json",
        [O_xs[0], O_ys[0]],
        [O_xs[2], O_ys[2]],
        O_xs,
        O_ys,
        O_lons,
        O_lats,
        proj_epsg3413,
    )

    speak.notquiet(args, "\nPlotting the grid.")
    plt.figure(1, figsize=(12, 10), dpi=150)
    plt.rc("text", usetex=True)
    plt.rc("font", family="serif")

    plt.subplot(1, 2, 1)  # EPSG:3413 in Bamber
    # NOTE: Basemap adds false Eastings and Northings dependent upon your corners.
    #      The proj4 projections, however, don't have any, so you'll always have
    #      to pass the map longs and lats.
    map_ll_lon, map_ll_lat = proj_eigen_gl04c(
        B_xs[0] - 800000.0, B_ys[0] - 1150000.0, inverse=True
    )
    map_ur_lon, map_ur_lat = proj_eigen_gl04c(
        B_xs[2] + 800000.0, B_ys[2] + 300000.0, inverse=True
    )
    glmap = Basemap(
        llcrnrlon=map_ll_lon,
        llcrnrlat=map_ll_lat,
        urcrnrlon=map_ur_lon,
        urcrnrlat=map_ur_lat,
        projection="stere",
        lat_ts=71.0,
        lat_0=90.0,
        lon_0=321.0,
        resolution="l",
    )

    glmap.fillcontinents(color="gray", lake_color="white")
    glmap.drawcoastlines()

    T_X, T_Y = glmap(T_lons, T_lats)
    glmap.plot(T_X, T_Y, "ko-", label="Bamber template grid")

    TE_X, TE_Y = glmap(TE_lons, TE_lats)
    glmap.plot(TE_X, TE_Y, "mo-", label="EPSG:3413 template grid")

    B_X, B_Y = glmap(B_lons, B_lats)
    glmap.plot(B_X, B_Y, "bo-", label="Bamber grid")

    E2B_X, E2B_Y = glmap(E_lons, E_lats)
    glmap.plot(E2B_X, E2B_Y, "ro-", label="EPSG:3413 grid")

    O2B_X, O2B_Y = glmap(O_lons, O_lats)
    glmap.plot(O2B_X, O2B_Y, "go-", label="Extended EPSG:3413 grid")

    plt.title("Greenland using the Bamber projection")
    plt.legend(loc="lower center", fancybox=True, shadow=True)

    plt.subplot(1, 2, 2)  # Bamber in EPSG:3413
    # NOTE: Basemap adds false Eastings and Northings dependent upon your corners.
    #      The proj4 projections, however, don't have any, so you'll alsways have
    #      to passs the map longs and lats.
    map_ll_lon, map_ll_lat = proj_epsg3413(
        E_xs[0] - 800000.0, E_ys[0] - 1150000.0, inverse=True
    )
    map_ur_lon, map_ur_lat = proj_epsg3413(
        E_xs[2] + 800000.0, E_ys[2] + 300000.0, inverse=True
    )
    glmap = Basemap(
        llcrnrlon=map_ll_lon,
        llcrnrlat=map_ll_lat,
        urcrnrlon=map_ur_lon,
        urcrnrlat=map_ur_lat,
        resolution="l",
        epsg=3413,
    )

    glmap.fillcontinents(color="gray", lake_color="white")
    glmap.drawcoastlines()

    T2E_X, T2E_Y = glmap(T_lons, T_lats)
    glmap.plot(T2E_X, T2E_Y, "ko-", label="Bamber template grid")

    TE_X, TE_Y = glmap(TE_lons, TE_lats)
    glmap.plot(TE_X, TE_Y, "mo-", label="EPSG:3413 template grid")

    B2E_X, B2E_Y = glmap(B_lons, B_lats)
    glmap.plot(B2E_X, B2E_Y, "bo-", label="Bamber grid")

    E_X, E_Y = glmap(E_lons, E_lats)
    glmap.plot(E_X, E_Y, "ro-", label="EPSG:3413 grid")

    O_X, O_Y = glmap(O_lons, O_lats)
    glmap.plot(O_X, O_Y, "go-", label="Extended EPSG:3413 grid")

    plt.title("Greenland using the EPSG:3413 projection")
    plt.legend(loc="lower center", fancybox=True, shadow=True)

    plt.tight_layout()
    plt.savefig("Bamber_v_EPSG3413.png", bbox_inches="tight")
    if args.show:
        plt.show()

    nc_base.close()
    speak.notquiet(args, "   Done!\n")


if __name__ == "__main__":
    main(parse_args())
