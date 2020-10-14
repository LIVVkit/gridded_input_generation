#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compare a newly generated CISM-like input dataset to an older version."""
import logging
import argparse
from pathlib import Path
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import seaborn as sns

__author__ = "Michael Kelleher"
plt.style.use("fivethirtyeight")


def common_test(ref, test, dtype=""):
    """Take two lists, compare contents."""
    common = set(ref).intersection(test)
    test_only = set(test).difference(ref)
    ref_only = set(ref).difference(test)
    log = logging.getLogger("test")

    if test_only:
        log.info(f"   {dtype} MISSING FROM  REF: {test_only}")
    if ref_only:
        log.info(f"   {dtype} MISSING FROM TEST: {ref_only}")

    return common, ref_only, test_only


def dict_diff(attrs_ref, attrs_test):
    """Find out what changed."""
    log = logging.getLogger("test")
    ref_vars = list(attrs_ref.keys())
    test_vars = list(attrs_test.keys())
    common, ref_only, test_only = common_test(ref_vars, test_vars, "Attribute")
    st_len = 90
    if common:
        for attr in common:
            _atr_ref = attrs_ref[attr]
            _atr_test = attrs_test[attr]
            if _atr_ref == _atr_test:
                if isinstance(_atr_ref, str) and len(_atr_ref) > st_len:
                    _atr_ref = f"{_atr_ref[:st_len]}..."

                log.info(f"   MATCH: {attr}: {_atr_ref}")
            else:
                # If they match is fine to shorten the attr, but if they're
                # different, print the whole thing
                log.info(
                    f"\n   ======== DIFFERENCE: {attr} =========="
                    f"\n      REF: {_atr_ref}"
                    f"\n     TEST: {_atr_test}\n"
                    "   ----------------------------------"
                )


def check_vars(ref, test):
    """Get a dictionary of metadata, matching reference and test files."""
    ref_vars = ref.data_vars
    test_vars = test.data_vars
    common, ref_only, test_only = common_test(ref_vars, test_vars, "Data Var")
    log = logging.getLogger("test")

    for dvar in sorted(common):
        log.info(f"\n{dvar.upper()} ATTRIBUTES")
        dict_diff(ref[dvar].attrs, test[dvar].attrs)

    return common, ref_only, test_only


def get_map_transform(dset, var):
    """Get appropriate map transform information."""
    try:
        grid_map = dset[var].grid_mapping
    except AttributeError:
        grid_map = None

    try:
        _map = dset[grid_map]
    except KeyError:
        grid_map = None
    if grid_map is not None:
        _map = dset[grid_map]
        params = {}
        _name = _map.attrs.get("grid_mapping_name", None)
        if _name == "polar_stereographic":
            proj = ccrs.Stereographic
            param_names = [
                (
                    "straight_vertical_longitude_from_pole",
                    "central_longitude",
                    0.0,
                ),
                ("latitude_of_projection_origin", "central_latitude", -90.0),
                ("standard_parallel", "true_scale_latitude", -71.0),
            ]
            for in_param, out_param, default in param_names:
                params[out_param] = _map.attrs.get(in_param, default)
        else:
            proj = ccrs.PlateCarree
        tform = proj(**params)
    elif (
        dset[var]["y1"].min() == -3400000
        and dset[var]["x1"].min() == -800000.0
    ):
        # Reference file for Bamber grid doesn't have grid map, and is on this
        # grid, so default to that if the lower left corner matches
        tform = ccrs.Stereographic(
            central_longitude=321.0,
            central_latitude=90.0,
            true_scale_latitude=71.0,
        )

    elif dset[var]["y1"].min() == -3333500.0 == dset[var]["x1"].min():
        # Reference file for Bamber grid doesn't have grid map, and is on this
        # grid, so default to that if the lower left corner matches
        tform = ccrs.Stereographic(
            central_longitude=0.0,
            central_latitude=-90.0,
            true_scale_latitude=71.0,
        )

    else:
        tform = ccrs.PlateCarree()

    return tform


def plot_single(data, varname, title, skip=1, axes=None):
    """Plot one contour panel, add boxplot if the axis input is None."""
    print(f"     {title}")
    tform = get_map_transform(data, varname)
    if axes is None:
        fig = plt.figure(figsize=(13, 7))
        axes = [
            fig.add_subplot(1, 2, 1, projection=tform),
            fig.add_subplot(1, 2, 2),
        ]
        single_var = True
    else:
        fig = None
        single_var = False
        axes = [axes]

    if (
        "missing_value" not in data[varname].attrs
        or "missing_value" not in data[varname].encoding
    ):
        _data = np.ma.masked_values(
            data[varname][0], data[varname][0].max()
        ).compressed()

    else:
        _data = data[varname][0]

    try:
        vmin, vmax = np.nanpercentile(_data, (2, 98))
    except TypeError:
        vmin, vmax = (-1, 1)

    if vmin < 0 < vmax:
        _allmax = np.max(np.abs([vmin, vmax]))
        vmin = -_allmax
        vmax = _allmax
        cmap = "RdBu_r"
    else:
        cmap = "viridis"

    _cf = axes[0].pcolormesh(
        data.x1[::skip],
        data.y1[::skip],
        data[varname][0, ::skip, ::skip],
        transform=tform,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )
    axes[0].set_title(title)
    plt.colorbar(_cf, ax=axes[0], pad=0.02, shrink=0.8)
    axes[0].gridlines()
    axes[0].coastlines(resolution="10m")

    if single_var:
        _data1d = np.ma.masked_invalid(
            data[varname].values.ravel()
        ).compressed()
        axes[1].boxplot(
            [_data1d], labels=[title], whis=[1, 99],
        )
        axes[1].text(
            1.1,
            np.nanpercentile(_data1d, 50),
            annote(data[varname].values.ravel()),
            horizontalalignment="left",
        )

        sns.despine(ax=axes[1], offset=10)
        fig.suptitle(varname, fontsize=18)

    return fig


def plot_sideby(ref, test, cmn_vars, out_dir, skip=1):
    """"Make side-by-side plot comparing reference to test."""

    for var in sorted(cmn_vars):
        if "epsg" in var or "mcb" in var or "mapping" in var:
            continue

        _diff = xr.Dataset({var: test[var] - ref[var]})
        _gridmap = test[var].attrs.get("grid_mapping")
        _diff[var] = _diff[var].assign_attrs({"grid_mapping": _gridmap})
        _diff[_gridmap] = test[_gridmap]

        _proj_lat = test[_gridmap].attrs.get(
            "latitude_of_projection_origin", 90
        )
        if _proj_lat > 0:
            proj = ccrs.NorthPolarStereo(central_longitude=0)
        else:
            proj = ccrs.SouthPolarStereo(central_longitude=0)

        fig = plt.figure(figsize=(13, 13))
        axes = [
            fig.add_subplot(2, 2, i + 1, projection=proj) for i in range(3)
        ]
        # No transform for the last axis
        axes.append(fig.add_subplot(2, 2, 4))

        print(f"Plot {var}")
        plot_single(ref, var, "Reference", skip, axes=axes[0])
        plot_single(test, var, "Test", skip, axes=axes[1])
        plot_single(_diff, var, "Test - Reference", skip, axes[2])

        print(f"     boxplot")
        # Because boxplot uses np.percentile to compute the box location,
        # We need to mask invalid data, and compress to just the non-masked
        # values, so that the percentiles don't come out as NaNs

        _masktest = np.ma.masked_greater(
            np.ma.masked_invalid(test[var].values.ravel()), 1e36
        ).compressed()

        _maskref = np.ma.masked_greater(
            np.ma.masked_invalid(ref[var].values.ravel()), 1e36
        ).compressed()

        axes[3].boxplot(
            [_maskref, _masktest], labels=["Reference", "Test"], whis=[1, 99],
        )
        axes[3].text(
            1.1,
            np.nanpercentile(_masktest, 50),
            annote(test[var].values.ravel()),
            horizontalalignment="left",
        )
        axes[3].text(
            2.1,
            np.nanpercentile(_maskref, 50),
            annote(ref[var].values.ravel()),
            horizontalalignment="left",
        )
        axes[3].set_title("Boxplot")
        print(" save figure")
        plt.suptitle(var, fontsize=18)
        fig.subplots_adjust(
            left=0.02,
            bottom=0.04,
            right=0.95,
            top=0.97,
            wspace=0.12,
            hspace=0.01,
        )
        sns.despine(ax=axes[3], offset=10)
        plt.savefig(Path(out_dir, f"plt_{var}_compare.{EXTN}"))
        plt.close()


def annote(data):
    """Get text annotation for boxplot, so masked values are evident."""
    str_out = ""

    if any(np.isnan(data)):
        str_out += "MASK NaN\n"
    if any(data > 1e36):
        str_out += "MASK >1e36"
    if str_out == "":
        str_out = "NO MASK"

    return str_out


def make_parse():
    # parse the command line arguments
    parser = argparse.ArgumentParser()  # -h or --help automatically included!

    parser.add_argument(
        "-i",
        "--island",
        type=str,
        default="ais",
        help="Island to compare (one of ais, gis)",
    )
    parser.add_argument(
        "-p",
        "--proj",
        type=str,
        default="epsg3031",
        help=(
            "Projection to compare (for island == ant: epsg3031, "
            "for island == gis: mcb or epsg3413"
        ),
    )

    return parser.parse_args()


def main():
    """Define and load two files, call plots."""
    args = make_parse()

    if args.island == "ais" and args.proj == "epsg3031":
        # Test Antarctica
        # ref_file = "ncs/antarctica_8km_2014_01_14.nc"
        # ref_file = "ncs/antarctica_1km_2017_05_03.nc"
        ref_file = "ncs/antarctica_8km_2020_03_19.nc"
        test_file = "ncs/antarctica_8km_2020_10_13.nc"
        out_dir = "ais_compare"
        log_out = "ais"

    elif args.island == "gis" and args.proj == "mcb":
        # Test Bamber grid
        ref_file = "complete/greenland_8km_2016_12_01.mcb.nc"
        # ref_file = "complete/greenland_8km_2020_03_04.mcb.nc"
        test_file = "complete/greenland_8km_2020_04_21.mcb.nc"
        out_dir = "bamber_compare"
        log_out = "gis_bamber"

    elif args.island == "gis" and args.proj == "epsg3413":
        # Test EPSG:3413 grid
        ref_file = "complete/greenland_8km_2017_06_27.epsg3413.nc"
        test_file = "complete/greenland_8km_2020_04_20.epsg3413.nc"
        out_dir = "epsg3413_compare"
        log_out = "gis_epsg"

    if not Path(out_dir).exists():
        Path(out_dir).mkdir()

    ref = xr.open_dataset(ref_file, decode_times=False).load()
    test = xr.open_dataset(test_file, decode_times=False).load()

    logging.basicConfig(
        level=logging.INFO,
        # format="%(asctime)s %(name)-12s %(levelname)-8s %(" "message)s",
        format="%(message)s",
        datefmt="%m-%d %H:%M:%S",
        filename=Path(out_dir, f"compare_{log_out}.log"),
        filemode="w",
    )

    log = logging.getLogger("test")
    log.info("-" * 30)
    log.info(f"Comparing\n  Test: {test_file}\n        to\n  Ref : {ref_file}")
    skip = 1
    cmn_vars, ref_vars, test_vars = check_vars(ref, test)
    plot_sideby(ref, test, cmn_vars, out_dir, skip=skip)

    for _refvar in ref_vars:
        if "epsg" in _refvar or "mcb" in _refvar or "mapping" in _refvar:
            continue
        print(f"Plot reference {_refvar}")
        _fig = plot_single(ref, _refvar, "Reference", skip=skip)
        _fig.savefig(Path(out_dir, f"plt_{_refvar}_ref_only.{EXTN}"))

    for _testvar in test_vars:
        if "epsg" in _testvar or "mcb" in _testvar or "mapping" in _testvar:
            continue
        print(f"Plot test {_testvar}")
        _fig = plot_single(test, _testvar, "Test", skip=skip)
        _fig.savefig(Path(out_dir, f"plt_{_testvar}_test_only.{EXTN}"))


EXTN = "png"
if __name__ == "__main__":
    main()
