#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compare a newly generated CISM-like input dataset to an older version."""
import logging
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs

matplotlib.use("Agg")
__author__ = "Michael Kelleher"


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
    ref_vars = attrs_ref.keys()
    test_vars = attrs_test.keys()
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
                if isinstance(_atr_ref, str) and len(_atr_ref) > st_len:
                    _atr_ref = f"{_atr_ref[:st_len]}..."
                if isinstance(_atr_test, str) and len(_atr_test) > st_len:
                    _atr_test = f"{_atr_test[:st_len]}..."
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
        log.info(f"\nCHECK {dvar} ATTRIBUTES")
        dict_diff(ref[dvar].attrs, test[dvar].attrs)

    return common, ref_only, test_only


def get_map_transform(dset, var):
    """Get appropriate map transform information."""
    try:
        grid_map = dset[var].grid_mapping
    except AttributeError:
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
    else:
        tform = ccrs.PlateCarree()

    return tform


def plot_single(data, varname, title, skip=1, axes=None):
    """Plot one contour panel, add boxplot if the axis input is None."""
    print(f"     {title}")
    tform = get_map_transform(data, varname)
    if axes is None:
        fig = plt.figure(figsize=(13, 13))
        axes = [
            fig.add_subplot(1, 2, 1, projection=tform),
            fig.add_subplot(1, 2, 2),
        ]
        single_var = True
    else:
        fig = None
        single_var = False
        axes = [axes]

    vmin, vmax = np.nanpercentile(data[varname][0], (2, 98))

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
        # levels=np.linspace(*np.nanpercentile(ref[var][0], (5, 95)), 15),
    )
    axes[0].set_title(title)
    plt.colorbar(_cf, ax=axes[0])

    if single_var:
        axes[1].boxplot(
            [np.ma.masked_invalid(data[varname].values.ravel()).compressed(),],
            labels=[title],
            whis=[1, 99],
        )
        fig.suptitle(varname)

    return fig


def plot_sideby(ref, test, cmn_vars, skip=1):
    """"Make side-by-side plot comparing reference to test."""
    proj = ccrs.SouthPolarStereo(central_longitude=0)

    for var in sorted(cmn_vars):
        if "epsg" in var:
            continue

        fig = plt.figure(figsize=(13, 13))
        axes = [
            fig.add_subplot(2, 2, i + 1, projection=proj) for i in range(3)
        ]
        # No transform for the last axis
        axes.append(fig.add_subplot(2, 2, 4))

        _diff = xr.Dataset({var: test[var] - ref[var]})
        _gridmap = test[var].attrs.get("grid_mapping")
        _diff[var] = _diff[var].assign_attrs({"grid_mapping": _gridmap})
        _diff[_gridmap] = test[_gridmap]

        print(f"Plot {var}")
        plot_single(ref, var, "Reference", skip, axes=axes[0])
        plot_single(test, var, "Test", skip, axes=axes[1])
        plot_single(_diff, var, "Test - Reference", skip, axes[2])

        for axis in axes[:-1]:
            axis.gridlines()
            axis.coastlines()

        print(f"     boxplot")
        # Because boxplot uses np.percentile to compute the box location,
        # We need to mask invalid data, and compress to just the non-masked
        # values, so that the percentiles don't come out as NaNs
        axes[3].boxplot(
            [
                np.ma.masked_invalid(test[var].values.ravel()).compressed(),
                np.ma.masked_invalid(ref[var].values.ravel()).compressed(),
            ],
            labels=["Test", "Ref"],
            whis=[1, 99],
        )
        axes[3].set_title("Boxplot")
        print(" save figure")
        plt.suptitle(var)
        plt.savefig(f"plt_{var}_compare.png")
        plt.close()


def main():
    """Define and load two files, call plots."""
    # ref_file = "/Users/25k/Data/ant/complete/antarctica_1km_2017_05_03.nc"
    ref_file = "ncs/antarctica_1km_2018_05_14.nc"
    test_file = "ncs/antarctica_1km_2020_02_10.nc"
    ref = xr.open_dataset(ref_file, decode_times=False).load()
    test = xr.open_dataset(test_file, decode_times=False).load()

    logging.basicConfig(
        level=logging.INFO,
        # format="%(asctime)s %(name)-12s %(levelname)-8s %(" "message)s",
        format="%(message)s",
        datefmt="%m-%d %H:%M:%S",
        filename="compare_1km_ais.log",
        filemode="w",
    )
    log = logging.getLogger("test")
    log.info("-" * 30)
    log.info(f"Comparing\n  Test: {test_file}\n        to\n  Ref : {ref_file}")
    skip = 1
    cmn_vars, ref_vars, test_vars = check_vars(ref, test)
    plot_sideby(ref, test, cmn_vars, skip=skip)
    for _refvar in ref_vars:
        print(f"Plot reference {_refvar}")
        _fig = plot_single(ref, _refvar, "Reference", skip=skip)
        _fig.savefig(f"plt_ref_{_refvar}.png")

    for _testvar in test_vars:
        print(f"Plot test {_testvar}")
        plot_single(test, _testvar, "Test", skip=skip)
        _fig.savefig(f"plt_test_{_testvar}.png")


if __name__ == "__main__":
    main()
