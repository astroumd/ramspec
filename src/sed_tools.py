import logging
import numpy as np
from astropy.io import fits

logger = logging.getLogger('dev')
logger.setLevel(logging.DEBUG)


def rect(y, x):
    """Integrate the first dimension of y over x.

    y can be multidimension arrays but only the first dimension is integrated.

    """

    dx = x[1:] - x[:-1]
    yleft = y[:-1, ...]
    ret = np.dot(np.transpose(yleft), dx)  # np.dot sums over the last
    # dimension
    return np.transpose(ret)


def from_f_lambda_to_f_nu(f_lambda, wavelengths, unit='micron'):
    """
    Use f_nu = f_lambda * lambda^2 / c to convert f_lambda to f_nu

    Args:
        f_lambda: (array) flux density
        wavelengths: (array)
        unit: (str) unit of the wavelengths, which should also be the unit
            of the 'per wavelength' in f_lambda

    Returns:
        nu: (array) frequency in Hz
        f_nu: (array) flux density in unit of ... per Hz

    """

    light_speed = c.c.to(f'{unit}/s').value  # in micron/s
    nu = light_speed / wavelengths      # Hz
    f_nu = f_lambda * wavelengths**2 / c
    return nu, f_nu


def from_f_nu_to_jansky(f_nu):
    """
    1 Jy = 10−23 erg⋅s−1⋅cm−2⋅Hz−1

    Args:
        f_nu: (array) f_nu in units of erg/s/cm2/Hz

    Returns:
        flux density in jansky
    """

    return 1e-23 * f_nu


def spectral_integrate(data, wavelength):
    """
    Integrate the spectral which is the first axis of the data cube.

    Returns:
        flux (array): the flux in units of erg/s/cm2/arcsec2. The 'per
            micron' or 'per Hz' is integrated and get away.
    """

    assert data.shape[0] == len(wavelength), \
        f"MyError: data.shape[0] = {data.shape[0]} but len(wavelength) = " \
        "{len(wavelength)}"
    return rect(data, wavelength)


def spacial_integrate(fn_fits, is_int_over_space=True):
    """Read data cube (where the first axis is wavelength or equivalent) from
    fits file and return spacially integrated SED.

    Args:
        fn_fits (str): filename of the fits data
        is_int (bool): default True. Toggl Integrate over solid angle and
        output in the unit W/m2/um or equivalent. Other wise output in unit
        W/m2/um/arcsec2 or equivalent

    Returns:
        _wavelengths (array): wavelength read from the fits
        _spectrum (array): the spacially integrated spectral energy density,
            Unit: if is_int, W/m2/um or equivalent; if not is_int: W/m2/um/arcsec2
            or equivalent.
    """

    with fits.open(fn_fits) as _hdus:
        _wavelengths = np.array([wavel[0] for wavel in _hdus[1].data])  # micron
        _data = _hdus[0].data
        _spectrum = np.mean(np.mean(_data, axis=-1), axis=-1)
        if is_int_over_space:
            hd = _hdus[0].header
            # assert hd['BUNIT'] == 'erg/s/cm2/arcsec2'
            solid_angle_per_pixel = hd['CDELT1'] * hd['CDELT2']  # arcsec2
            logger.info("Units:" + str(hd['BUNIT']+'*'+hd['CUNIT1']+'*'+hd['CUNIT2']))
            _spectrum *= solid_angle_per_pixel * (_data.shape[1] * _data.shape[2])
            print("Spactiall integrated asusming a source distance shown in "
                  "the respective .ski file")
        else:
            logger.info("Units:" + str(_hdus[0].header['BUNIT']))
    return _wavelengths, _spectrum

def get_sed(prefix, kind="sed", is_int=0):
    """OLD, DO NOT USE!
    Read SED from fits or .dat file.
    Given a fits data with grids of spectra, return a spatially integrated spectrum

    Args
        prefix (str): prefix to get the data filename prefix_total.fits or prefix_sed.dat
        kind (str): default "sed". Either "sed" or "fits"
        is_int (bool): default 0. Effective when kind = "fits". Integrate over solid
            angle and output in the unit W/m2/um. Other wise output in unit W/m2/um/arcsec2.

    Returns:
        _wavelengths (array): wavelength read from the fits
        _spectrum (array): the spacially integrated spectral energy density,
        in unit of W/m2/um if is_int else W/m2/um/arcsec2

    """

    assert kind in ["sed", "fits"]
    if kind == "fits":
        fn = f'{prefix}_total.fits'
        return spacial_integrate(fn, is_int_over_space=is_int)
    elif kind == "sed":
        fn = f'{prefix}_sed.dat'
        _wavelengths, _spectrum = np.loadtxt(fn, unpack = True)
        with open(fn, 'r') as _f:
            _f.readline()
            logger.info("Units:" + _f.readline().split(' ')[-1])
    return _wavelengths, _spectrum


def combine_spec(specs, wavelengths, line_str, line_wavelength, line_width='resol'):
    """Combine pixels of spectrum with pixels of single-line surface brightness.

    Note that the units of spec2 / wavelength2 should be the same as that of
    specs.

    Parameters
    ----------
    specs: (m, n1, n2, ...) array
        The pixeled spectrum
    wavelengths: (m, ) array
        List of wavelengths
    spec2: (n1, n2, ...) array
        The single-line surface brightness (e.g. Ly-alpha, H-alpha). Should have the same units
        as specs EXPECT lacking a wavelength component (W/m2/arcsec2)
    wavelength2 (double)
        The wavelength of spec2 (same units as wavelengths)
    width2 ('resol' or double):
        (Default: 'resol') The width of spec2. If set to resol, use the
        resolution at the corresponding wavelength in wavelengths.

    """

    assert specs.shape[1:] == line_str.shape, "specs: {}, spec2: {}".format(
        specs.shape, line_str.shape)
    pick = np.argmax(wavelengths >= line_wavelength)
    assert pick > 0 and pick < len(wavelengths) - 1, "line_wavelength is not in the range of wavelengths"
    pick -= 1
    dlam_instru = wavelengths[pick + 1] - wavelengths[pick]  # micron
    ret = specs.copy()
    if isinstance(line_width, str):
        if line_width == "resol":
            logger.warn('Warning: Using spectral resolution as the width of H-alpha. This is not physical unless the spectral resolution matches real instrument.')
            ret[pick, ...] += line_str / dlam_instru
        else:
            raise SystemExit("width2 should be either 'resol' or a float.")
    else:
        pick_c = np.argmax(wavelengths >= line_wavelength)
        dl = wavelengths[pick_c] - wavelengths[pick_c - 1]
        if line_width < dl:
            ret[pick_c, ...] += line_str / line_width
        else:
            lam_l, lam_r = line_wavelength - line_width/2, line_wavelength + line_width/2
            pick_l = np.argmax(wavelengths >= lam_l)
            pick_r = np.argmax(wavelengths >= lam_r)
            assert pick_l and pick_r  # both of them should not be 0
            if pick_r == pick_l:
                pick_r += 1
            ret[pick_l:pick_r, ...] += line_str / line_width
    return ret


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
        # print(f"Applying filter with l_pick = {l_pick}, return data[pick, ...] / width_at_l_pick * l_pick")
        return data[pick, ...] / (wavelengths[pick] - wavelengths[pick - 1]) * l_pick
    if l_c is not None:
        l_min = l_c - l_w / 2
        l_max = l_c + l_w / 2
    left = np.argmax(wavelengths >= l_min)
    right = np.argmax(wavelengths >= l_max)
    if right == left:
        right += 1
    ret = np.zeros(data.shape[1:])
    for i in range(left, right):
        ret += data[i, ...] * (wavelengths[i] - wavelengths[i - 1])
    # print(f"Applying filter with l_min = {l_min}, l_max = {l_max}, width ="
    #       f" {l_max - l_min}, # of bins: {right-left}")
    # print(f"Applying filter The actual wavelength range is from {wavelengths[left]} to "\
    #       f"{wavelengths[right]}, width ="
    #       f" {wavelengths[right] - wavelengths[left]}, # of bins: {right-left}")
    if scale is not None:
        ret *= scale
    return ret
