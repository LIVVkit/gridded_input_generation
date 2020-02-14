# encoding: utf-8


import os

import numpy as np
import pandas as pd


def antarctica(data_location):
    drainage_numbers = list(range(1, 28))
    drainage = pd.read_csv(
        os.path.join(
            data_location, "Ant_Grounded_DrainageSystem_Polygons.txt"
        ),
        sep="\s+",
        names=["lat", "lon", "id"],
        header=None,
        skiprows=7,
    )

    icesheet_names = {28: "AP", 29: "WAIS", 30: "EAIS"}
    icesheet = pd.read_csv(
        os.path.join(data_location, "Ant_Grounded_IceSheets_Polygon-1.txt"),
        sep="\s+",
        names=["lat", "lon", "id"],
        header=None,
        skiprows=9,
    )

    for n in drainage_numbers:
        basin = drainage.loc[drainage[id] == n, ["lat", "lon"]]
