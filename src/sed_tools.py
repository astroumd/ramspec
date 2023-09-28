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


def spacial_integrate(fn_fits, is_int=True):
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
            in unit of W/m2/um or equivalent if is_int else W/m2/um/arcsec2
            or equivalent.
    """

    with fits.open(fn_fits) as _hdus:
        _wavelengths = np.array([wavel[0] for wavel in _hdus[1].data])  # micron
        _data = _hdus[0].data
        _spectrum = np.mean(np.mean(_data, axis=-1), axis=-1)
        if is_int:
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
        return spacial_integrate(fn, is_int=is_int)
    elif kind == "sed":
        fn = f'{prefix}_sed.dat'
        _wavelengths, _spectrum = np.loadtxt(fn, unpack = True)
        with open(fn, 'r') as _f:
            _f.readline()
            logger.info("Units:" + _f.readline().split(' ')[-1])
    return _wavelengths, _spectrum


def combine_spec(specs, wavelengths, spec2, wavelength2, width2='resol'):
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
        (Default: None) The width of spec2. If set to None (default), use the
        resolution at the corresponding wavelength in wavelengths.

    """

    assert specs.shape[1:] == spec2.shape, "specs: {}, spec2: {}".format(
        specs.shape, spec2.shape)
    pick = np.argmax(wavelengths >= wavelength2)
    dlam_instru = wavelengths[pick + 1] - wavelengths[pick]  # micron
    ret = specs.copy()
    if width2 is "resol":
        logger.warn('Warning: Using resol of the wavelengths. This is not physical unless '\
                       'wavelengths represents real instrument.')
        ret[pick, ...] += spec2 / dlam_instru
    else:
        assert type(width2) is float, "width2 should be a float but is {}".format(type(width2))
        # if width2 > dlam_instru:
        #     raise SystemExit("Your width2 is larger than the width of the instrument."\
        #                      "This is not supported right now.")
        lam_l, lam_r = wavelength2 - width2/2, wavelength2 + width2/2
        pick_l = np.argmax(wavelengths >= lam_l)
        pick_r = np.argmax(wavelengths >= lam_r)
        assert pick_l and pick_r  # either of them should not be 0
        if pick_r == pick_l:
            pick_r += 1
        ret[pick_l:pick_r, ...] += spec2 / width2
    return ret

