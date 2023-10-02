import os
import numpy as np
import matplotlib.pyplot as plt
import json
import scienceplots

plt.style.use(['science', 'no-latex'])


def plot_halpha(data_dir, outs, plot_dir, box_fraction=1.0):
    """Make mock photo.
    
    args:
        data_dir (str): path to the data directory
        outs (list of integers): list of output numbers
        plot_dir (str): path to the plot directory
        feature (str): prefix of the output files from SKIRT simulation.
    """

    ha_data_dir = os.path.join(data_dir, "halpha")

    radius = box_fraction / 2

    for out in outs:
        # H-alpha, unit: erg s-1 cm-2 arcsec-2
        surfb_ext = np.load(f"{ha_data_dir}/sb_with_dust_out{out:05d}.npy")
        surfb_no_dust = np.load(f"{ha_data_dir}/sb_no_dust_out{out:05d}.npy")

        f, ax = plt.subplots(figsize=(6, 6))
        extent = [-radius, radius, -radius, radius]
        im = ax.imshow(np.log10(surfb_ext)[::-1, :], extent=extent,
                        vmin=-18, vmax=-12, cmap='Reds')
        ax.set(xlabel='y (code unit)', ylabel='z (code unit)', )
        ax.set_title('H-alpha flux per sr')
        cb = f.colorbar(im, ax=ax)
        cb.set_label(r'log F (erg/s/cm2/sr)')
        plt.savefig(f"{plot_dir}/halpha-with-dust-out{out:05d}.png")

        f, ax = plt.subplots(figsize=(6, 6))
        extent = [-radius, radius, -radius, radius]
        im = ax.imshow(np.log10(surfb_no_dust)[::-1, :], extent=extent,
                        vmin=-18, vmax=-12, cmap='Reds')
        ax.set(xlabel='y (code unit)', ylabel='z (code unit)', )
        ax.set_title('H-alpha flux (without dust) per sr')
        cb = f.colorbar(im, ax=ax)
        cb.set_label(r'log F (erg/s/cm2/sr)')
        plt.savefig(f"{plot_dir}/halpha-no-dust-out{out:05d}.png")


def plot_combined_spec(data_dir, outs, plot_dir, feature):
    """Plot combined spectrum.

    Args:
        data_dir (str): path to the data directory
        outs (list of integers): list of output numbers
        feature (str): prefix of the output files from SKIRT simulation.
    """
    
    combined_data_dir = os.path.join(data_dir, "combined", feature)

    for out in outs:
        # fo_combined = f"{combined_data_dir}/combined_spec_out{out:05d}.json"
        # with open(fo_combined, 'r') as fin:
        #     combined = json.load(fin)
        # wavel = np.array(combined['wavelength'])
        # spec = np.array(combined['spec'])

        data = np.loadtxt(f"{combined_data_dir}/combined_spectral_density_with_dust_out{out:05d}.txt", skiprows=2)
        wavel = data[:, 0]
        spec = data[:, 1]

        f, ax = plt.subplots()
        plt.plot(wavel, spec)
        ax.set(xscale='log', xlabel=r"$\lambda$ ($\mu$m)", yscale='log',
                ylabel=r'$F_{\lambda}$ [erg/(s cm$^2$ arcsec$^2$ $\mu$m)]',
                ylim=[6e-16, 6e-11],
        )
        plt.savefig(f"{plot_dir}/sed_with_ha_out{out:05d}.pdf")
        print(f"{plot_dir}/sed_with_ha_out{out:05d}.pdf saved")
