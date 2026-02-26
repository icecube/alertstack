import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u


def prepare_ras_aitoff(ras):
    """Adjust right ascensions for aitoff projection

    Parameters
    ----------
    ras: `numpy.array[float]`
        right ascensions (in radiants)
    """
    ras_aitoff = np.empty(len(ras))
    ras_aitoff[ras<=np.pi] = ras[ras<=np.pi]
    ras_aitoff[ras>np.pi] = ras[ras>np.pi] - 2*np.pi
    return ras_aitoff

def prepare_gal_coords(gal_coords):
    """Get some galactic coordinates and return the
    correspondent right ascension and declination.

    Parameters
    ----------
    gal_coords: `numpy.array[astropy.coordinates.SkyCoord]`
        The galactic coordinates.
    """
    ra_vals = gal_coords.icrs.ra.rad
    dec_vals = gal_coords.icrs.dec.rad
    ra_vals_aitoff = prepare_ras_aitoff(ra_vals)
    return ra_vals_aitoff, dec_vals

def plot_gp_coords(gp_cut):
    """Given a cut in galactic latitude, plot this cut.

    Parameters
    ----------
    gp_cut: `float`
        Cut in galactic latitude.
    """
    if gp_cut > 0.:
        gp_lons = np.linspace(0., 360., 100)
        gal_p = SkyCoord(
            l = gp_lons*u.deg, b = gp_cut*u.deg, frame='galactic'
        )
        gal_m = SkyCoord(
            l = gp_lons*u.deg, b = -gp_cut*u.deg, frame='galactic'
        )
        ra_vals_aitoff_p, dec_vals_p = prepare_gal_coords(gal_p)
        ra_vals_aitoff_m, dec_vals_m = prepare_gal_coords(gal_m)
        plt.plot(
            1,
            1,
            color="black",
            linestyle='dotted',
            label=f"|gal lat| = {int(gp_cut)} deg",
            linewidth=2
        )
        plt.scatter(ra_vals_aitoff_p, dec_vals_p, s=1.5, c="black")
        plt.scatter(ra_vals_aitoff_m, dec_vals_m, s=1.5, c="black")

def plot_dec_cut(dec_cut):
    """Given a cut in declination, plot this cut.

    Parameters
    ----------
    dec_cut: `float`
        Cut in declination.
    """
    if dec_cut > -90.:
        ras_bord = np.linspace(-np.pi, np.pi, 100)
        des_bord = np.empty(len(ras_bord))
        des_bord.fill(dec_cut*np.pi/180.)
        plt.plot(
            ras_bord, 
            des_bord, 
            color="black", 
            label=f"Dec. = {int(dec_cut)} deg", 
            linewidth=2
        )

def get_weight_name(weighter):
    """Given an hypothesis, extract the name for the weight

    Parameters
    ----------
    weighter: `alertstack.Hypothesis`
        The hypothesis.
    """
    weight_name = ""
    weight_name_components = weighter.name.split("_")
    for i, word in enumerate(weight_name_components):
        weight_name += word[0].upper() + word[1:]
        if i + 1 < len(weight_name_components):
            weight_name += " "
    if weighter.unit is not None:
        weight_name += f" [{weighter.unit}]"
    return weight_name

def plot_pdf(catalogue_obj, nside=16, final_pdf=False):
    """Count the number of sources per pixel and plot the result.

    Parameters
    ----------
    nside: `int`
        nside of the binning
    final_pdf: `bool`
        If True, use the final pdf used for that catalog.
        If False, just do the binning.
    """
    hd_nside = 256
    if final_pdf:
        nside = 128
        bins_probs = hp.ud_grade(catalogue_obj.bkg_distribution, hd_nside)
    else:
        catalogue = catalogue_obj.parse_data()
        bins_per_source = hp.ang2pix(
            nside,
            np.pi/2. - catalogue["dec_rad"],
            catalogue["ra_rad"]
        )
        bins = np.arange(hp.nside2npix(nside)+1)
        counts_per_bin, _ = np.histogram(bins_per_source, bins=bins)
        bins_probs = counts_per_bin/np.sum(counts_per_bin)
    bins_probs = bins_probs / hp.nside2pixarea(nside)
    bins_probs = hp.ud_grade(bins_probs, hd_nside)
    bins_cothetas, bins_phis = hp.pix2ang(
        hd_nside, np.arange(hp.nside2npix(hd_nside))
    )
    bins_thetas = np.pi / 2. - bins_cothetas
    phis_aitoff = prepare_ras_aitoff(bins_phis)
    im = plt.scatter(
        phis_aitoff, bins_thetas, c=bins_probs, cmap="Blues", s=0.01
    )
    cb = plt.gcf().colorbar(im)
    cb.ax.tick_params(labelsize="large")
    cb.set_label(label="Density of Probability [rad-2]", fontsize="large")

def plot_catalogue(
    catalogue_obj,
    s,
    scatter_label,
    title,
    weighter=None,
    tw=False,
    nside=None,
    final_pdf=False,
    scramble=False,
    scramble_size=20,
):
    """Plot a catalogue with an aitoff projection

    Parameters
    ----------
    catalogue_obj: `alertstack.ScrambleCatalogue`
        Object for the catalogue to plot
    s: `int`
        Size of dots in the scatter plot
    scatter_label: `str`
        Label for the sources (to show in the legend)
    title: `str`
        Title for the plot
    weighter: `alertstack.Hypothesis | None`
        Hypothesis to weight the catalogue. If None, ignore weights
    tw: `bool`
        If the weight requires or not a time window
    nside: `int | None`
        Plot binning of sources with the desired nside
    scramble: `bool`
        Show the example of a scramble
    """

    catalogue = catalogue_obj.parse_data()

    plt.figure(figsize=(10,4))
    ax = plt.subplot(111, projection='aitoff')

    ras_aitoff = prepare_ras_aitoff(catalogue["ra_rad"])

    if nside is not None or final_pdf:
        plot_pdf(catalogue_obj, nside=nside, final_pdf=final_pdf)

    if weighter is not None and tw:
        weights = weighter.weight_catalogue(
            catalogue, nu_at=0., ignore_times=True
        )
    elif weighter is not None:
        weights = weighter.weight_catalogue(catalogue)
    else:
        weights = "tab:orange"
        cmap = None
    im = plt.scatter(
        ras_aitoff,
        catalogue["dec_rad"],
        s=s,
        label=scatter_label,
        c=weights,
        norm="log",
        cmap="plasma"
    )
    if scramble:
        scrambled_cat = catalogue_obj.scramble()
        plt.scatter(
            prepare_ras_aitoff(scrambled_cat["ra_rad"]),
            scrambled_cat["dec_rad"],
            c="tab:blue",
            s=scramble_size,
            label="Scramble",
            marker="x",
        )
    plot_gp_coords(catalogue_obj.gp_threshold)
    plot_dec_cut(catalogue_obj.min_declination)

    plt.legend(loc=(-0.15, 0.95), fontsize='large')
    if weighter is not None:
        cb = plt.gcf().colorbar(im)
        cb.ax.tick_params(labelsize="large")
        cb.set_label(label=get_weight_name(weighter), fontsize="large")

    plt.grid()

    plt.xticks(fontsize="large")
    plt.yticks(fontsize="large")

    plt.title(title)
    plt.xlabel("R.A. [deg]", fontsize="large")
    plt.ylabel("Dec. [deg]", fontsize="large")