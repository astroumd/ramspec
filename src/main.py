import os
import numpy as np
import argparse
import json
import cv2
from astropy.io import fits
from sed_tools import combine_spec, get_sed
from ramtools import ramses, utilities, to_skirt
from ramtools.cacherunbase import CacheRun
from h_alpha import halpha_sb
from ramspec.src.amr2cube import exe_amr2cube
# from plotit import abmag
import sed_tools


DATA = "../data/new_2022"
SAM = "/startrek/chongchong/Projects/2017-RAMSES"
FILTER = os.path.join(os.path.dirname(__file__), "filters")
LAM_HALPHA = 0.65628  # micron


def combine_fits_with_ha(hdul, ha, fo=None, R=None, use_spectral_width=False):
    """ hdul: fits data of the stellar continuum. The dimensions of the data are (wavelength, x,
    y)

    Args:
        hdul (fits data): continuum data, should have the unit
            erg/s/cm2/micron/arcsec2 (OLD: W/m2/micron/arcsec2)
        ha (numpy array): h-alpha, should have the unit
            erg/s/cm2/arcsec2 (OLD: W/m2/arcsec2)
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
        data[pick, ...] = R * ha    # erg/s/cm2/arcsec2
    if fo is not None:
        hdul.writeto(fo, overwrite=True)
    # spec2 = ha               # W/m2/arcsec2
    # return combine_spec(data, wavelengths, spec2, lam_halpha, width2='resol')
    return data


def image_filter(data, wavelengths, l_c=None, l_w=None, l_min=None,
                 l_max=None, l_pick=None, scale=None):
    """ Given a data cube (3-D numpy array), return the filtered image in 2-D

    Args:
        data: (3-D numpy array) spectrum in pixels, unit:
            erg/s/cm2/micron/arcsec2 or equivalent
        wavelengths: (numpy array) wavelengths, the first dimension of data
        l_c: (float) the center of the filter
        l_w: (float) the width of the filter
        l_min: (float) the left boundary of the filter
        l_max: (float) the right boudnary of the filter
        l_pick: (float) if not None, overwrite all previous four parameters.
            Pick a single image at the given wavelength and return data *
            l_pick.
        scale: (float) scale up or down the spectral energy

    Returns:
        (2-D numpy array)
    """

    if l_pick is not None:
        pick = np.argmax(wavelengths >= l_pick)
        return data[pick, ...] * l_pick
    if l_c is not None:
        l_min = l_c - l_w / 2
        l_max = l_c + l_w / 2
    left = np.argmax(wavelengths >= l_min)
    right = np.argmax(wavelengths >= l_max)
    ret = np.zeros(data.shape[1:])
    for i in range(left, right):
        ret += data[i, ...] * (wavelengths[i+1] - wavelengths[i])
    print(f"The actual wavelength range is from {wavelengths[left]} to "\
          f"{wavelengths[right]}, width ="
          f" {wavelengths[right] - wavelengths[left]}, # of bins: {right-left}")
    if scale is not None:
        ret *= scale
    return ret


def run_sub_region(args):

    # # 'crab-small'
    # center = (0.32061501, 0.49432314, 0.45865159)
    # width = 0.063035        # box width in boxlen unit, crab_small
    # jobid = "2.2.2.v2"
    # out = 19
    # axis = 1
    # # .dat files are ready to use
    # halpha_fns = "/startrek2nb/chongchong/Sam/test_amr2cube/max_{f}_l13.dat"    # crab_small
    # skirt = "sam19_max_small_1e9ph_total_total.fits"  # crab-small

    # 'crab'
    # center = (0.32061501, 0.49432314, 0.45865159)
    # width = 0.1008558850  # box width in boxlen unit, crab
    # jobid = "2.2.2.v2"
    # out = 19
    # axis = 2
    # lma = 12
    # is_redo_amr2cube = False
    # skirt = "sam19_max_1e9ph_sedsb_total.fits"  # crab, axis z

    jobdir = args.jobdir
    out = args.out
    # diri = f"{SAM}/Job{jobid}/output_{out:05d}"
    diri = f"{jobdir}/output_{out:05d}"
    center = [.5, .5, .5] if args.center is None else args.center
    width = args.width
    axis = args.axis
    lma = args.lma
    skirt = args.skirt

    #-------------------------------------------------
    # Prepare cube data using amr2cube for H-alpha calculation.
    #-------------------------------------------------
    print(">>>>>>>>>> Preparing cube data using amr2cube...")
    halpha_fns = "amr2cube/{field}.dat"
    os.makedirs("amr2cube", exist_ok=True)
    # CacheRun(exe_amr2cube, os.path.abspath(__file__))(
    #     diri, center, width, lma, fo_fmt=halpha_fns)
    exe_amr2cube(diri, center, width, lma, fo_fmt=halpha_fns, is_redo=args.redo_amr2cube)
    print("<<<<<<<<< Done")
    # r = ramses.Ramses(jobid=jobid, ram_dir=SAM)

    #-------------------------------------------------
    # prepare particle data
    #-------------------------------------------------
    dir_skirt = "skirt"
    os.makedirs(dir_skirt, exist_ok=1)
    to_skirt.to_skirt(None, out, f"{dir_skirt}/part_CK", family="CK",
                      center=center, width_boxlen=width, letdie=False,
                      skip_exist=True, jobdir=jobdir)
    if (not args.do_ha_only) and (skirt is None):
        return

    #-------------------------------------------------
    # Compute H-alpha
    #-------------------------------------------------
    print(">>>>>>>>>> Computing H-alpha...")
    boxlen = utilities.read_quant_from_ramses_info(
        f"{diri}/info_{out:05d}.txt", "boxlen"
    )
    unit_l = utilities.read_quant_from_ramses_info(
        f"{diri}/info_{out:05d}.txt", "unit_l"
    )
    # width = width * r.unit_l   # cm
    width_cm = width * boxlen * unit_l   # cm
    den = halpha_fns.format(field='den')
    xHII = halpha_fns.format(field='xHII')
    surfb = halpha_sb(den, xHII, width_cm, axis=axis)  # erg s-1 cm-2 arcsec-2
    # data_ha = np.loadtxt(fn_ha)     # erg s-1 cm-2 arcsec-2
    sb_SI_to_cgs = 1000.            # from W/m2/arcsec2 to erg/s/cm2/arcsec2
    # data_ha_SI = data_ha / sb_SI_to_cgs      # W/m2/arcsec2
    # data_ha_SI = surfb / sb_SI_to_cgs      # W/m2/arcsec2
    # save ha data into fits file
    print("<<<<<<<<< Done\n")
    if args.do_ha_only:
        dir_ha = "ha"
        os.makedirs(dir_ha, exist_ok=True)
        fits.PrimaryHDU(surfb).writeto("ha/ha.fits", overwrite=True)
        return
    if skirt is None:
        return

    #-------------------------------------------------
    # Final step: combine h-alpha with dust continuum
    #-------------------------------------------------
    print(">>>>>>>>>> Combining h-alpha with dust continuum...")
    # Combine h-alpha with dust continuum
    # fn_ha = f"{DATA}/{thedir}/ha/{ha}"
    diro = "skirt+ha"
    os.makedirs(diro, exist_ok=True)
    skirt_base = os.path.basename(skirt)
    fn_dust_and_ha = f"{diro}/{skirt_base.replace('.fits', '+ha.fits')}"
    fn_dust_and_ha_intensity = fn_dust_and_ha.replace('.fits', '-indensity.fits')
    # fn_dust = f"skirt/{skirt}"
    with fits.open(skirt) as hdul:
        data = hdul[0].data
        wavelengths = np.array(
            [wavel[0] for wavel in hdul[1].data])  # micron
        # 1. convert hdul unit from W/m2/arcsec2/micron to erg/s/cm2/arcsec2
        # by multiplying with lambda. 2. Merge data_ha * R with hdul. R
        # is the resolution parameter (R = lambda / d_lambda)
        R = args.R     # spectra resolution = lambda / d_lambda
        # combine_fits_with_ha(hdul, data_ha_SI, fo=fn_dust_and_ha, R=R)
        # rescale surfb to make it match SKIRT fits
        s_ha = surfb.shape[0]   # first dimension
        s_skirt = data.shape[1]   # first dimension
        if s_ha != s_skirt:
            if (int(round(s_skirt / s_ha)) == 2) and \
                    args.use_repeat_as_resize:
                # simply repeat
                t = np.repeat(np.repeat(surfb, 2, axis=0), 2, axis=1)
                surfb = np.zeros(data.shape[1:])
                surfb[:t.shape[0], :t.shape[1]] = t
            else:
                surfb = cv2.resize(surfb, data.shape[1:],
                                   interpolation=cv2.INTER_LINEAR)
        combine_fits_with_ha(hdul, surfb, fo=fn_dust_and_ha, R=R)
        # hdul.writeto(fn_dust_and_ha, overwrite=True)

    with fits.open(skirt) as hdul:
        # write intensity (micron-1) to fits
        data_comb = combine_fits_with_ha(
            hdul, surfb, fo=fn_dust_and_ha_intensity, use_spectral_width=True)

    # write out r g b fits data
    # fil_g = {"l_c": 6455/1e4, "l_w": 82/1e4}    # in micron
    # fil_b = {"l_c": 5013/1e4, "l_w": 47/1e4}
    # fil_g = {"l_c": 6455/1e4, "l_w": 200/1e4}    # in micron
    # fil_g = {"l_min": 5500/1e4, "l_max": 6500/1e4}    # in micron
    # # fil_b = {"l_c": 5013/1e4, "l_w": 300/1e4}
    # fil_b = {"l_c": 1500/1e4, "l_w": 70/1e4}
    fn_filter = f"{FILTER}/{args.filter}.json"
    with open(fn_filter, 'r') as f:
        filter = json.load(f)
        fil_r = "ha" if "fil_r" not in filter.keys() else filter[
            "fil_r"]
        fil_g = filter['fil_g']
        fil_b = filter['fil_b']
    dir_rgb = f"rgb-{args.filter}"
    os.makedirs(dir_rgb, exist_ok=1)
    fnr = f"{dir_rgb}/red.fits"
    fng = f"{dir_rgb}/green.fits"
    fnb = f"{dir_rgb}/blue.fits"
    # from W/m2/micron/arcsec2 to erg/s/micron/cm2/arcsec2
    sb_SI_to_cgs = 1000.
    dataset = {}
    for fil, color in zip([fil_r, fil_g, fil_b], ['r', 'g', 'b']):
        if isinstance(fil, str):
            if fil == "ha":
                dataset[color] = surfb * R
        else:
            dataset[color] = image_filter(
                data_comb, wavelengths, **fil) * sb_SI_to_cgs
    # da_r = surfb * R
    # da_g = image_filter(data, wavelengths, **fil_g) * sb_SI_to_cgs
    # da_b = image_filter(data, wavelengths, **fil_b) * sb_SI_to_cgs
    # print("red:", da_r.min(), da_r.max(), np.log10(da_r).mean())
    # print("green:", da_g.min(), da_g.max(), np.log10(da_g).mean())
    # print("blue:", da_b.min(), da_b.max(), np.log10(da_b).mean())
    # unit: erg/s/cm2/arcsec2
    fits.PrimaryHDU(dataset['r']).writeto(fnr, overwrite=True)
    fits.PrimaryHDU(dataset['g']).writeto(fng, overwrite=True)
    fits.PrimaryHDU(dataset['b']).writeto(fnb, overwrite=True)

    # # calculate AB magnitude
    # # I abandoned this because this does not make sence. The resulting
    # # numbers are too small because it is 'per arcsec2'
    # abm = abmag(data_comb, wavelengths)
    # fo = f"{diro}/ABmag.fits"
    # fits.PrimaryHDU(abm).writeto(fo, overwrite=True)

    # calculate spectrum-integrated map
    integrated = sed_tools.spectral_integrate(data_comb, wavelengths)
    hdu = fits.PrimaryHDU(integrated)
    hdu.header['UNIT'] = "erg/s/cm2/arcsec2"
    fo = f"{diro}/integrated-map.fits"
    hdu.writeto(fo, overwrite=True)

    # integrate over wavelength to get SED
    l, fl_times_l = sed_tools.spacial_integrate(fn_dust_and_ha)
    # now, fl_times_l = f_lambda * lambda, has unit 'erg/s/cm2'
    arr = np.array([l, fl_times_l]).T
    fo_sed = f"{diro}/SED.txt"
    np.savetxt(fo_sed, arr, header="unit1: micron\nunit2: erg/s/cm2")

    # # pick single wavelength
    # l_g = 6455/1e4
    # # l_b = 5013/1e4
    # l_b = 1500/1e4
    # dir_rgb = "rgb-single-bin"
    # os.makedirs(dir_rgb, exist_ok=1)
    # fnr = f"{dir_rgb}/red.fits"
    # fng = f"{dir_rgb}/green.fits"
    # fnb = f"{dir_rgb}/blue.fits"
    # data = hdul[0].data
    # wavelengths = np.array(
    #     [wavel[0] for wavel in hdul[1].data])  # micron
    # # from W/m2/micron/arcsec2 to erg/s/micron/cm2/arcsec2
    # sb_SI_to_cgs = 1000.
    # da_r = surfb * R
    # da_g = image_filter(data, wavelengths, l_pick=l_g) * \
    #        sb_SI_to_cgs
    # da_b = image_filter(data, wavelengths, l_pick=l_b) * \
    #        sb_SI_to_cgs
    # for name, d_i in zip(['red', 'green', 'blue'],
    #                      [da_r, da_g, da_b]):
    #     print(name + ':', d_i.min(), d_i.max(), np.percentile(
    #         d_i, [.9, .99, 0.99999]))
    # # unit: erg/s/cm2/arcsec2
    # fits.PrimaryHDU(da_r).writeto(fnr, overwrite=True)
    # fits.PrimaryHDU(da_g).writeto(fng, overwrite=True)
    # fits.PrimaryHDU(da_b).writeto(fnb, overwrite=True)

        # # try converting axis3 to energy (eV)
    # with fits.open(fn_dust_and_ha) as hdul:
    #     print("[0] header")
    #     print(hdul[0].header)
    #     print("[1] header")
    #     print(hdul[1].header)
    #     hdul[0].header['CUNIT3'] = 'eV'
    #     hdul[1].header['TUNIT1'] = 'eV'
    #     print("[0] header")
    #     print(hdul[0].header)
    #     print("[1] header")
    #     print(hdul[1].header)
    #     hdul.writeto(f"{diro}/try.fits", overwrite=1)

    print("<<<<<<<<< Done")
    return


def compute_sed_for_a_snap(jobid, out, lmax, width, axis, R, letdie=True):

    # jobid = "2.2.2"
    # out = 20
    # lmax = 9
    # width = 0.8
    # axis = 1
    # R = 1000
    skirt = f"../data/yorp07/run_v6/Job{jobid}/out{out}/main_CK_corrected_nd_1e6ph_sedsb_total.fits"
    if letdie:
        skirt = skirt.replace("_nd_", "_")
    halpha_fns = f"../data/data_amr2cube/Job{jobid}/out{out}_{{" \
                 f"field}}_l{lmax}.dat"
    strdie = "letdie-" if letdie else ""
    fo_sed = f"../data/local_run_v6/Job{jobid}/sed-{strdie}R{R}/out{out}.txt"
    fo_dust_and_ha = f"../data/local_run_v6/Job{jobid}/skirt_and_ha-{strdie}" \
                     f"R{R}/out{out}.fits"
    if not os.path.exists(skirt):
        print(f"Skipping {jobid}, out{out} because the following file does not "
              f"exist:", skirt)
        return
    if os.path.exists(fo_sed):
        print(f"Skipping because {fo_sed} exists")
        return
    os.makedirs(os.path.dirname(fo_sed), exist_ok=1)
    os.makedirs(os.path.dirname(fo_dust_and_ha), exist_ok=1)

    #-------------------------------------------------
    # H-alpha
    #-------------------------------------------------
    diri = f"{SAM}/Job{jobid}/output_{out:05d}"
    boxlen = utilities.read_quant_from_ramses_info(
        f"{diri}/info_{out:05d}.txt", "boxlen"
    )
    unit_l = utilities.read_quant_from_ramses_info(
        f"{diri}/info_{out:05d}.txt", "unit_l"
    )
    width_cm = width * boxlen * unit_l   # cm
    den = halpha_fns.format(field='den')
    xHII = halpha_fns.format(field='xHII')
    surfb = halpha_sb(den, xHII, width_cm, axis=axis)  # erg s-1 cm-2 arcsec-2

    #-------------------------------------------------
    # Combine stellar continuum with H-alpha
    #-------------------------------------------------
    with fits.open(skirt) as hdul:
        combine_fits_with_ha(hdul, surfb, fo=fo_dust_and_ha, R=R)
    l, fl_times_l = sed_tools.spacial_integrate(fo_dust_and_ha)
    # now, fl_times_l = f_lambda * lambda, has unit 'erg/s/cm2'
    arr = np.array([l, fl_times_l]).T
    np.savetxt(fo_sed, arr, header="unit1: micron\nunit2: erg/s/cm2")
    return


def compute_sed_for_snaps(jobid, snaps, **kwargs):
    for snap in snaps:
        compute_sed_for_a_snap(jobid, snap, **kwargs)


def compute_sed_for_all(letdie=True):
    kwargs = dict(
        lmax = 9,
        width = 0.8,
        axis = 1,
        R = 1000,
        letdie = letdie,
    )
    compute_sed_for_snaps('2.2.2', range(16, 49 + 1, 2), **kwargs)
    compute_sed_for_snaps('3.2.2', range(14, 44 + 1, 2), **kwargs)
    compute_sed_for_snaps('4.2.1', range(14, 48 + 1, 2), **kwargs)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare SKIRT simualtions and make H-alpha data")
    parser.add_argument("jobdir", nargs='?',
                        help="dirctory of the RAMSES job")
    parser.add_argument("out", nargs='?', type=int,
                        help="dirctory of the RAMSES job")
    parser.add_argument("width", nargs='?', type=float,
                        help="width in boxlen unit")
    parser.add_argument("lma", nargs='?', type=int,
                        help="maximum level of refinement for amr2cube")
    parser.add_argument("-a", "-axis", dest="axis", type=int, default=1,
                        help="projection axis")
    parser.add_argument("-c", "-center", dest='center', nargs='+', type=float,
                        help="(list of 3) the center in boxlen unit")
    parser.add_argument("-redo-amr2cube", action="store_true",
                        help="toggle force redo amr2cube")
    parser.add_argument("-s", "-skirt-fits", dest="skirt",
                        help="(optional) the SKIRT fits data")
    parser.add_argument("-do-ha-only", action="store_true",
                        help="toggle compute and save H-alpha only")
    parser.add_argument("-R", "-spectra-resolution", dest="R",
                        type=float, default=100,
                        help="spectra resolution (= lambda / d_lambda)")
    parser.add_argument("-f", "-filter", dest="filter",
                        help="path to the filter file")
    parser.add_argument("-use-repeat-as-resize", action="store_true",
                        help="Toggle use repeating as image resizing")
    parser.add_argument("-task", default="region",
                        help="Specify the task to conduct. Default is None "
                             "in which case main() is executed.")
    parser.add_argument("-letdie", action="store_true",
                        help="toggle let stars die (default: False)")
    return parser.parse_args()


if __name__ == "__main__":

    args = parse_args()
    if args.task == "region":
        exit(run_sub_region(args))
    elif args.task == 'sed':
        exit(compute_sed_for_all(args.letdie))