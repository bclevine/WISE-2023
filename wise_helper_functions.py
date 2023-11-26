import numpy as np
import pandas as pd
from astropy.coordinates import Distance
import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib import cm

iso = pd.read_csv("Gaia_iso.csv")


def get_isochrone(age, metallicity, distance):
    # Loads in an isochrone from a file given age (in logspace), metallicity, and distance (in astropy units)
    z = f"p{metallicity:.2f}" if metallicity >= 0 else f"m{np.abs(metallicity):.2f}"
    distmod = Distance(distance).distmod.value
    g = iso[f"G_{age:.2f}_{z}"] + distmod
    bp = iso[f"BP_{age:.2f}_{z}"] + distmod
    rp = iso[f"RP_{age:.2f}_{z}"] + distmod
    b = iso[f"B_{age:.2f}_{z}"] + distmod
    v = iso[f"V_{age:.2f}_{z}"] + distmod
    r = iso[f"R_{age:.2f}_{z}"] + distmod
    return (
        g[~np.isnan(g)],
        bp[~np.isnan(bp)],
        rp[~np.isnan(rp)],
        b[~np.isnan(b)],
        v[~np.isnan(v)],
        r[~np.isnan(r)],
    )


def plot_positions(catalog):
    color = "G" if "G" in catalog.columns else "R"
    plt.figure(figsize=(10, 7))
    plt.scatter(catalog["ra"], catalog["dec"], c=catalog[color], cmap="viridis_r")
    plt.xlabel("Right Ascension [degrees]")
    plt.ylabel("Declination [degrees]")
    plt.colorbar(label=f"{color} [magnitude]")
    plt.show()


def plot_CMD(catalog, opacity=0.5, **kwargs):
    mag = "G" if "G" in catalog.columns else "R"
    col = "B-V" if "B-V" in catalog.columns else "BP-RP"
    plt.figure(figsize=(10, 7))
    plt.scatter(catalog[col], catalog[mag], alpha=opacity, label="Stars", **kwargs)
    plt.gca().invert_yaxis()
    plt.xlabel(f"{col} [magnitude]")
    plt.ylabel(f"{mag} [magnitude]")
