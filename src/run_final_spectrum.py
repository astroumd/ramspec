#!/startrek/chongchong/anaconda3/envs/yt/bin/python -u

""" Make final spectrum, combining continuum with H-alpha """

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# import astropy.units as u, astropy.constants as c
import json
import f90nml
from astropy.io import fits

from academicpython import plottools as pt
from . import h_alpha
from .sed_tools import get_sed, combine_spec

matplotlib.rcParams['figure.dpi'] = 300
plt.style.use(['science', 'no-latex'])

C_hu_halpha = 3.026e-12  # in erg, = 1.89 eV
C_arcsec2_to_sr = 2.35044e-11
lam_halpha = 0.65628  # um


def make_combined_spec(data_dir, outs, feature):
    """Make combined spectrum, combining continuum with H-alpha.

    Args:
        data_dir (str): path to the data directory
        outs (list of integers): list of output numbers
        feature (str): prefix of the output files from SKIRT simulation.
    """
    
    # load halpha data from file
    skirt_data_dir = os.path.join(data_dir, "skirt")
    ha_data_dir = os.path.join(data_dir, "halpha")
    combined_data_dir = os.path.join(data_dir, "combined", feature)

    os.makedirs(combined_data_dir, exist_ok=True)

    for out in outs:
        # H-alpha
        data_json = json.load(open(os.path.join(ha_data_dir, f"line_str_out{out:05d}.json")))
        ha_strength = data_json['halpha_strength']
        ha_unit = data_json['unit']
        assert(ha_unit == "erg s-1 cm-2 arcsec-2")
        # surfb_ext = np.load(f"{ha_data_dir}/sb_with_dust_out{out:05d}.npy")
        # surfb_no_dust = np.load(f"{ha_data_dir}/sb_no_dust_out{out:05d}.npy")


        # fn_fits = f"{skirt_data_dir}/{feature}_total.fits"
        # with fits.open(fn_fits) as hdul:
        #     combine_fits_with_ha(hdul, surfb, fo=fo_dust_and_ha, R=R)

        #-------- spatially integrated spectrum of continuum --------#
        root_path = f"{skirt_data_dir}/out{out:05d}"
        fn_fits = f"{root_path}/{feature}_total.fits"
        wavel, spec = get_sed(f'{root_path}/{feature}')
        wavel2, spec2 = get_sed(f'{root_path}/{feature}', kind="fits", is_int=1)

        # integrate over space
        with fits.open(fn_fits) as _hdus:
            hd = _hdus[0].header
            solid_angle_per_pixel = hd['CDELT1'] * hd['CDELT2']  # arcsec2
            total_arcsec2 = hd['NAXIS1'] * hd['NAXIS2'] * solid_angle_per_pixel
            assert hd['CUNIT1'] == "arcsec" and hd['CUNIT2'] == "arcsec"
        spec_per_arcsec2 = spec / total_arcsec2   # W / (m2 arcsec2 micron)
        spec_per_arcsec2 *= 1000.0                # erg / (s cm2 arcsec2 micron)
        # dic['dust_wavelength'][out] = NoIndent(list(wavel))
        # dic['dust_spec'][out] = NoIndent(list(spec_per_arcsec2))

        # Combine H-alpha with dust extinction
        lam_halpha = 0.65628   # um
        ha_width = lam_halpha / 100  # l/dl = 100
        halpha_sb_mean = ha_strength 
        halpha_sb_mean_arr = np.array([halpha_sb_mean]).reshape([1,1])  # make it 2D, 1x1
        spec_per_arcsec2_mat = spec_per_arcsec2[:, np.newaxis, np.newaxis] # make it 3D, Nx1x1
        data_with_ha = combine_spec(
            spec_per_arcsec2_mat, wavel, halpha_sb_mean_arr,
            lam_halpha, width2=ha_width)
        
        combined = {
            'wavelength': wavel,
            'wavelength_unit': 'micron',
            'spec': data_with_ha[:, 0, 0],
            'spec_unit': 'erg*/*(s*cm2*arcsec2*micron)',
        }
        # convert ndarray in combined to list
        for k, v in combined.items():
            if isinstance(v, np.ndarray):
                combined[k] = v.tolist()
        # save combined spectrum
        fo_combined = f"{combined_data_dir}/combined_spec_out{out:05d}.json"
        with open(fo_combined, 'w') as fout:
            json.dump(combined, fout, indent=2)

            
def plot_combined_spec(data_dir, outs, plot_dir, feature):
    """Plot combined spectrum.

    Args:
        data_dir (str): path to the data directory
        outs (list of integers): list of output numbers
        feature (str): prefix of the output files from SKIRT simulation.
    """
    
    combined_data_dir = os.path.join(data_dir, "combined", feature)

    for out in outs:
        fo_combined = f"{combined_data_dir}/combined_spec_out{out:05d}.json"
        with open(fo_combined, 'r') as fin:
            combined = json.load(fin)
        wavel = np.array(combined['wavelength'])
        spec = np.array(combined['spec'])
        f, ax = plt.subplots()
        plt.plot(wavel, spec)
        ax.set(xscale='log', xlabel=r"$\lambda$ ($\mu$m)", yscale='log',
                ylabel=r'$F_{\lambda}$ [erg/(s cm$^2$ arcsec$^2$ $\mu$m)]',
                ylim=[6e-16, 6e-11],
        )
        plt.savefig(f"{plot_dir}/sed_with_ha_out{out:05d}.pdf")
        print(f"{plot_dir}/sed_with_ha_out{out:05d}.pdf saved")

