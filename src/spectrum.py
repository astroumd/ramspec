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
import cv2
import shutil

from academicpython import plottools as pt
from . import h_alpha
from .sed_tools import get_sed, combine_spec, image_filter, spacial_integrate

THISDIR = os.path.dirname(os.path.realpath(__file__))

matplotlib.rcParams['figure.dpi'] = 300
# plt.style.use(['science', 'no-latex'])

C_hu_halpha = 3.026e-12  # in erg, = 1.89 eV
C_arcsec2_to_sr = 2.35044e-11
LAM_HALPHA = 0.65628  # um


def make_combined_spec(data_dir, outs, feature, filter="filter1", R=100):
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
        ha_with_dust_strength = data_json['halpha_with_dust_strength_mean']
        assert(data_json['unit'] == "erg s-1 cm-2 arcsec-2")
        surfb_ext = np.load(f"{ha_data_dir}/sb_with_dust_out{out:05d}.npy")
        surfb_no_dust = np.load(f"{ha_data_dir}/sb_no_dust_out{out:05d}.npy")

        # fn_fits = f"{skirt_data_dir}/{feature}_total.fits"
        # with fits.open(fn_fits) as hdul:
        #     combine_fits_with_ha(hdul, surfb, fo=fo_dust_and_ha, R=R)

        #-------- spatially integrated spectrum of continuum --------#
        root_path = f"{skirt_data_dir}/out{out:05d}"
        fn_fits = f"{root_path}/{feature}_total.fits"
        wavel, spec = get_sed(f'{root_path}/{feature}')
        wavel2, spec2 = get_sed(f'{root_path}/{feature}', kind="fits", is_int=1)

        def is_close(arr1, arr2):
            fails = []
            for i in range(len(arr1)):
                if arr1[i] is np.nan:
                    if arr2[i] is np.nan:
                        continue
                    else:
                        fails.append(i)
                        continue
                if np.abs(arr1[i]) < 1e-20 and np.abs(arr2[i]) < 1e-20:
                    continue
                if np.abs(spec[i] / spec2[i]) > 3 or np.abs(spec[i] / spec2[i]) < 1/3:
                    fails.append(i)
            return fails

        # check the convergence of spec and spec2. They have to be close.
        assert abs(spec.max() / spec2.max() - 1.0) < 0.01
        fails_ = is_close(spec, spec2)
        if sum(fails_) > 0:
            raise ValueError("spec and spec2 are not close at {fails_}")

        # integrate over space
        with fits.open(fn_fits) as _hdus:
            hd = _hdus[0].header
            solid_angle_per_pixel = hd['CDELT1'] * hd['CDELT2']  # arcsec2
            total_arcsec2 = hd['NAXIS1'] * hd['NAXIS2'] * solid_angle_per_pixel
            assert hd['CUNIT1'] == "arcsec" and hd['CUNIT2'] == "arcsec"
            assert _hdus[0].header['BUNIT'] == 'W/m2/micron/arcsec2'

        spec_per_arcsec2 = spec / total_arcsec2   # W / (m2 arcsec2 micron)
        spec_per_arcsec2 *= 1000.0                # erg / (s cm2 arcsec2 micron)

        # Combine H-alpha with dust extinction and get spectral density, F_lambda, and save to file
        ha_width = LAM_HALPHA / R  # l/dl = R
        halpha_sb_mean = ha_with_dust_strength 
        halpha_sb_mean_arr = np.array([halpha_sb_mean]).reshape([1,1])  # make it 2D, 1x1
        spec_per_arcsec2_mat = spec_per_arcsec2[:, np.newaxis, np.newaxis] # make it 3D, Nx1x1
        data_with_ha = combine_spec(
            spec_per_arcsec2_mat, wavel, line_str=halpha_sb_mean_arr,
            line_wavelength=LAM_HALPHA, line_width=ha_width)
        arr = np.array([wavel, data_with_ha[:, 0, 0]]).T
        np.savetxt(
            f"{combined_data_dir}/combined_spectral_density_with_dust_out{out:05d}.txt",
            arr, 
            header="unit1: micron\nunit2: erg/(s*cm2*micron*arcsec2)")
        
        # write combined data cube to fits
        with fits.open(fn_fits) as _hdus:
            data = _hdus[0].data
            s_ha = surfb_ext.shape[0]   # first dimension
            s_skirt = data.shape[1]   # first dimension
            if s_ha != s_skirt:
                surfb_ext = cv2.resize(surfb_ext, data.shape[1:], interpolation=cv2.INTER_LINEAR)
            combine_fits_with_ha_new(_hdus, surfb_ext, fo=f"{combined_data_dir}/combined_data_cube_with_dust_out{out:05d}.fits", R=R)
        with fits.open(fn_fits) as _hdus:
            data = _hdus[0].data
            s_ha = surfb_no_dust.shape[0]   # first dimension
            s_skirt = data.shape[1]   # first dimension
            if s_ha != s_skirt:
                surfb_no_dust = cv2.resize(surfb_no_dust, data.shape[1:], interpolation=cv2.INTER_LINEAR)
            combine_fits_with_ha_new(_hdus, surfb_no_dust, fo=f"{combined_data_dir}/combined_data_cube_no_dust_out{out:05d}.fits", R=R)
        # spec_density = f"{combined_data_dir}/combined_data_cube_with_dust_spectral_density_out{out:05d}.fits"
        # with fits.open(fn_fits) as _hdus:
        #     data = _hdus[0].data
        #     s_ha = surfb_ext.shape[0]   # first dimension
        #     s_skirt = data.shape[1]   # first dimension
        #     if s_ha != s_skirt:
        #         surfb_ext = cv2.resize(surfb_ext, data.shape[1:], interpolation=cv2.INTER_LINEAR)
        #     combine_fits_with_ha(_hdus, surfb_ext, fo=spec_density, use_spectral_width=True)

        
def apply_filters(data_dir, outs, feature, R=100, filter="filter1"):

    combined_data_dir = os.path.join(data_dir, "combined", feature)
    out_dir = os.path.join(data_dir, "filtered-image", feature)
    os.makedirs(out_dir, exist_ok=1)

    for out in outs:

        with fits.open(f"{combined_data_dir}/combined_data_cube_with_dust_out{out:05d}.fits") as hdul:
            surfb = hdul[0].data
            wavelengths = np.array([wavel[0] for wavel in hdul[1].data])  # micron
        
        # write out r g b fits data
        # fil_g = {"l_c": 6455/1e4, "l_w": 82/1e4}    # in micron
        # fil_b = {"l_c": 5013/1e4, "l_w": 47/1e4}
        # fil_g = {"l_c": 6455/1e4, "l_w": 200/1e4}    # in micron
        # fil_g = {"l_min": 5500/1e4, "l_max": 6500/1e4}    # in micron
        # # fil_b = {"l_c": 5013/1e4, "l_w": 300/1e4}
        # fil_b = {"l_c": 1500/1e4, "l_w": 70/1e4}
        fn_filter = f"{THISDIR}/../tools/filters/{filter}.json"
        with open(fn_filter, 'r') as f:
            filter_dict = json.load(f)
        color_names = {'r': 'red', 'g': 'green', 'b': 'blue'}
        for fil in filter_dict.keys():
            params = filter_dict[fil]
            if isinstance(params, str):
                if params == "ha":
                    params = {"l_pick": LAM_HALPHA}
            ds = image_filter(surfb, wavelengths, **params) # * sb_SI_to_cgs
            fits.PrimaryHDU(ds).writeto(f"{out_dir}/out{out:05d}_{color_names[fil]}.fits", overwrite=True)
        
        # copy filter file to out_dir and overwrite
        shutil.copy(fn_filter, os.path.join(out_dir, "filter.json"))

        # # # calculate AB magnitude
        # # # I abandoned this because this does not make sence. The resulting
        # # # numbers are too small because it is 'per arcsec2'
        # # abm = abmag(data_comb, wavelengths)
        # # fo = f"{diro}/ABmag.fits"
        # # fits.PrimaryHDU(abm).writeto(fo, overwrite=True)

        # # calculate spectrum-integrated map
        # integrated = sed_tools.spectral_integrate(data_comb, wavelengths)
        # hdu = fits.PrimaryHDU(integrated)
        # hdu.header['UNIT'] = "erg/s/cm2/arcsec2"
        # fo = f"{diro}/integrated-map.fits"
        # hdu.writeto(fo, overwrite=True)


def combine_fits_with_ha(hdul, ha, fo=None, R=None, use_spectral_width=False):
    """ Combine fits data with H-alpha. Will modify the data in hdul. Will change the unit of the data from W/m2/micron/arcsec2 to erg/s/cm2/micron/arcsec2.

    Args:
        hdul (fits data): continuum data, should have the unit W/m2/micron/arcsec2.
            The dimensions of the data are (wavelength, x, y)
        ha (numpy array): h-alpha, should have the unit
            erg/s/cm2/arcsec2 (OLD: W/m2/arcsec2). The integrated indensity, F_l * dl
        fo (str): If not None, will write into fits file with the given
            filename.
        R (float): spectrum resolution
        use_spectral_width (bool): toggle use spectral width as the
            resolution and return intensity (F_lambda) instead of F

    Return:
        None
    """

    assert hdul[0].header['CUNIT3'] == 'micron'
    assert hdul[0].header['BUNIT'] == 'W/m2/micron/arcsec2', \
        f"BUNIT = {hdul[0].header['BUNIT']}"
    # assert hdul[0].header['BUNIT'] == "erg/s/cm2/micron/arcsec2", f"BUNIT = {hdul[0].header['BUNIT']}"
    wavelengths = np.array([wavel[0] for wavel in hdul[1].data])  # micron
    data = hdul[0].data     # W/m2/micron/arcsec2
    sb_SI_to_cgs = 1000.  # from W/m2/arcsec2 to erg/s/cm2/arcsec2
    data *= sb_SI_to_cgs    # erg/s/cm2/micron/arcsec2
    hdul[0].header['BUNIT'] = "erg/s/cm2/micron/arcsec2"
    pick = np.argmax(wavelengths >= LAM_HALPHA)
    if use_spectral_width:
        data[pick, ...] += ha / (wavelengths[pick+1] - wavelengths[pick])
    else:   # do I = I * lambda
        for i in range(len(wavelengths)):
            data[i, ...] *= wavelengths[i]     # erg/s/cm2/arcsec2
        hdul[0].header['BUNIT'] = "erg/s/cm2/arcsec2"
        if ha.shape[0] != data.shape[1]:
            if ha.shape[0] - data.shape[1] == 1:
                ha = ha[:-1, :-1]
            elif ha.shape[0] - data.shape[1] == -1:
                data = data[:, :-1, :-1]
            else:
                raise ValueError(f"MyError: ha and skirt data have number of pixels different by more than 1. ha.shape = {ha.shape}, data.shape = {data.shape}")
        # data[pick, ...] = R * ha    # erg/s/cm2/arcsec2
        # the old one has the wrong unit. 
        data[pick, ...] = ha / (LAM_HALPHA / R)    # erg/s/cm2/micron/arcsec2
    if fo is not None:
        hdul.writeto(fo, overwrite=True)
    # spec2 = ha               # W/m2/arcsec2
    # return combine_spec(data, wavelengths, spec2, lam_halpha, width2='resol')
    # return data


def combine_fits_with_ha_new(hdul, ha, fo, R):
    """ Combine fits data with H-alpha. Will modify the data in hdul. Will change the unit of the data from W/m2/micron/arcsec2 to erg/s/cm2/micron/arcsec2.

    Args:
        hdul (fits data): continuum data, should have the unit W/m2/micron/arcsec2.
            The dimensions of the data are (wavelength, x, y)
        ha (numpy array): h-alpha, should have the unit
            erg/s/cm2/arcsec2 (OLD: W/m2/arcsec2). The integrated indensity, F_l * dl
        fo (str): If not None, will write into fits file with the given
            filename.
        R (float): spectrum resolution
        use_spectral_width (bool): toggle use spectral width as the
            resolution and return intensity (F_lambda) instead of F

    Return:
        None
    """

    assert hdul[0].header['CUNIT3'] == 'micron'
    assert hdul[0].header['BUNIT'] == 'W/m2/micron/arcsec2', \
        f"BUNIT = {hdul[0].header['BUNIT']}"
    wavelengths = np.array([wavel[0] for wavel in hdul[1].data])  # micron
    data = hdul[0].data     # W/m2/micron/arcsec2
    assert data.ndim == 3, f"data.ndim = {data.ndim}"
    assert ha.ndim == 2, f"ha.ndim = {ha.ndim}"
    assert data.shape[1] == ha.shape[0], f"data.shape = {data.shape}, ha.shape = {ha.shape}"
    sb_SI_to_cgs = 1000.  # from W/m2/arcsec2 to erg/s/cm2/arcsec2
    data *= sb_SI_to_cgs    # erg/s/cm2/micron/arcsec2
    hdul[0].header['BUNIT'] = "erg/s/cm2/micron/arcsec2"

    combined = combine_spec(data, wavelengths, line_str=ha, line_wavelength=LAM_HALPHA, line_width=LAM_HALPHA/R)

    data[...] = combined[...]
    if fo is not None:
        hdul.writeto(fo, overwrite=True)

    return

    pick = np.argmax(wavelengths >= LAM_HALPHA)
    if use_spectral_width:
        data[pick, ...] += ha / (wavelengths[pick+1] - wavelengths[pick])
    else:   # do I = I * lambda
        for i in range(len(wavelengths)):
            data[i, ...] *= wavelengths[i]     # erg/s/cm2/arcsec2
        hdul[0].header['BUNIT'] = "erg/s/cm2/arcsec2"
        if ha.shape[0] != data.shape[1]:
            if ha.shape[0] - data.shape[1] == 1:
                ha = ha[:-1, :-1]
            elif ha.shape[0] - data.shape[1] == -1:
                data = data[:, :-1, :-1]
            else:
                raise ValueError(f"MyError: ha and skirt data have number of pixels different by more than 1. ha.shape = {ha.shape}, data.shape = {data.shape}")
        # data[pick, ...] = R * ha    # erg/s/cm2/arcsec2
        # the old one has the wrong unit. 
        data[pick, ...] = ha / (LAM_HALPHA / R)    # erg/s/cm2/micron/arcsec2
    if fo is not None:
        hdul.writeto(fo, overwrite=True)
    # spec2 = ha               # W/m2/arcsec2
    # return combine_spec(data, wavelengths, spec2, lam_halpha, width2='resol')
    # return data

