#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Compare a newly generated CISM-like input dataset to an older version."""
import logging
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy.crs as ccrs
import matplotlib

matplotlib.use("Agg")
__author__ = "Michael Kelleher"


def common_test(ref, test, dtype=""):
    """Take two lists, compare contents."""
    common = set(ref).intersection(test)
    test_only = set(test).difference(ref)
    ref_only = set(ref).difference(test)
    log = logging.getLogger("test")

    if test_only:
        log.info(f"{dtype} MISSING FROM REF: {test_only}")
    if ref_only:
        log.info(f"{dtype} MISSING FROM TEST: {ref_only}")

    return common, ref_only, test_only


def dict_diff(attrs_ref, attrs_test):
    """Find out what changed."""
    log = logging.getLogger("test")
    ref_vars = attrs_ref.keys()
    test_vars = attrs_test.keys()
    common, ref_only, test_only = common_test(ref_vars, test_vars, "Attribute")
    if common:
        for attr in common:
            if attrs_ref[attr] == attrs_test[attr]:
                log.info(f"   MATCH: {attr}: {attrs_ref[attr]}")
            else:
                log.info(
                    f"\n   ======== DIFFERENCE: {attr} =========="
                    f"\n     REF: {attrs_ref[attr]}"
                    f"\n     TEST: {attrs_test[attr]}\n"
                    "   ----------------------------------"
                )


def check_vars(ref, test):
    """Get a dictionary of metadata, matching reference and test files."""
    ref_vars = ref.data_vars
    test_vars = test.data_vars
    common, ref_only, test_only = common_test(ref_vars, test_vars, "Data Var")
    log = logging.getLogger("test")

    for dvar in common:
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


def plot_sideby(ref, test, cmn_vars):
    """"Make side-by-side plot comparing reference to test."""
    proj = ccrs.SouthPolarStereo(central_longitude=0)

    for var in sorted(cmn_vars):
        if "epsg" in var:
            continue

        fig = plt.figure(figsize=(13, 5))
        axes = [
            fig.add_subplot(1, 2, i + 1, projection=proj) for i in range(2)
        ]
        tform = get_map_transform(ref, var)
        print(f"Plot {var} Ref")
        vmin, vmax = np.nanpercentile(ref[var][0], (5, 95))
        cf_ref = axes[0].pcolormesh(
            ref.x1,
            ref.y1,
            ref[var][0],
            transform=tform,
            vmin=vmin,
            vmax=vmax,
            # levels=np.linspace(*np.nanpercentile(ref[var][0], (5, 95)), 15),
        )
        plt.colorbar(cf_ref, ax=axes[0])
        print(f"Plot {var} Test")
        vmin, vmax = np.nanpercentile(test[var][0], (5, 95))
        cf_test = axes[1].pcolormesh(
            test.x1,
            test.y1,
            test[var][0],
            transform=tform,
            vmin=vmin,
            vmax=vmax
            # levels=np.linspace(*np.nanpercentile(test[var][0], (5, 95)), 15),
        )
        plt.colorbar(cf_test, ax=axes[1])
        plt.suptitle(var)
        for axis in axes:
            axis.gridlines()
            axis.coastlines()
        print("  save figure")
        plt.savefig(f"plt_{var}_compare.png")
        plt.close()


def main():
    """Define and load two files, call plots."""
    ref_file = "ncs/antarctica_1km_2018_05_14.nc"
    test_file = "ncs/antarctica_1km_2020_01_27.nc"
    ref = xr.open_dataset(ref_file, decode_times=False)
    test = xr.open_dataset(test_file, decode_times=False)

    logging.basicConfig(
        level=logging.INFO,
        # format="%(asctime)s %(name)-12s %(levelname)-8s %(" "message)s",
        format="%(message)s",
        datefmt="%m-%d %H:%M:%S",
        filename="compare_1km_ais.log",
        filemode="w",
    )

    cmn_vars, ref_vars, test_vars = check_vars(ref, test)
    plot_sideby(ref, test, cmn_vars)


if __name__ == "__main__":
    main()
