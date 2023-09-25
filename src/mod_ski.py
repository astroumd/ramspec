#!/startrek/chongchong/anaconda3/envs/skirt/bin/python -u

###!/usr/bin/env python
"""mod_ski.py
Modify .ski file and (may) run skirt using the python package pts.

"""

import os, sys
import logging
import matplotlib as mpl
# sys.path.append("/startrek/chongchong/PTS")
import pts.simulation as sm
# import pts.utils as ut
import pts.visual as vis
import pts.do
import f90nml
from astropy import units
from ramtools.utilities import read_quant_from_ramses_info

# sys.path.append("/n/startrek2nb/chongchong/Sam/coding/pkgs/utils")
# from tools import read_quant_from_ramses_info

pts.do.initializePTS()
mpl.use('Agg')
mpl.rcParams['figure.dpi'] = 300

def run_ski(ski, feature="main_CK"):
    """
    Warning: abandoned! Do NOT use unless rewrite this.
    nph (float): number of photons
    lmax (int or str): l_max used to determine the spacial resoltuion of the
        instrument
    fn_par (str): particle file name
    fov (float): field of view, in pc
    fn_ski (str): the ski file based on which we create the new ski file
    feature (str): name of the output .ski file. The final ski filename will be
        {feature}_#ph.ski where # is n_photon.

    """

    is_run = is_mod
    if not is_mod and not os.path.isfile(ski + "_parameters.xml"):
        is_run = True
    print("is_run =", is_run)
    if is_run:
        # perform the simulation
        logging.info("Executing " + skifn)
        skirt = sm.Skirt(path="/n/startrek/chongchong/SKIRT2/release/SKIRT/main/skirt")
        simulation = skirt.execute(ski, console='regular',)
    else:
        # if the simulation has already been performed, use this instead
        simulation = sm.createSimulation(prefix=ski)

    # plot the SED
    micron = sm.unit("micron")
    vis.plotSeds(simulation, #minWavelength=0.1*micron, maxWavelength=1000*micron,
                 decades=4, figSize=(7, 4.5), outFileName=f"{ski}_sed.png", )


def mod_ski(fi, fo, nph, lmax, fn_par, fn_nml, fn_ramses_info, version="2020"):
    """
    fi: input ski file
    fo: output ski file
    nph (float): number of photons
    lmax (int or str): l_max used to determine the spacial resoltuion of the
        instrument in order to make it consistant with my H-alpha program
    fn_par (str): particle file name
    fn_nml: RAMSKI namelist file
    fov (float): field of view, in pc
    """

    lmax = int(lmax)
    skifn = fo
    if os.path.isfile(skifn):
        if os.path.getmtime(skifn) > os.path.getmtime(fi):
            print("{fo} exists and is newer than this python script. Not doing anything")
            return 1
    ski = sm.SkiFile(fi)
    ski.setNumPrimaryPackets(nph)
    ski.setStringAttribute(
        "MonteCarloSimulation/sourceSystem/SourceSystem/sources/ParticleSource",
        "filename", fn_par)
    # ski.setStringAttribute(
    #     "MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument",
    #     "instrumentName", "sedsb")
    # set SED to CK
    # did by hand
    # calculate ny and ny
    nml = f90nml.read(fn_nml)
    side = float(nml["PARAMS"]["xmax"]) - float(nml["PARAMS"]["xmin"])
    sidey = float(nml["PARAMS"]["ymax"]) - float(nml["PARAMS"]["ymin"])
    sidez = float(nml["PARAMS"]["zmax"]) - float(nml["PARAMS"]["zmin"])
    assert sidey == side and sidez == side, \
            "Current version of mod_ski.py only works for equal-side box"
    nx = int(side * 2**lmax)
    ny = nx
    # fov = float(fov)
    # ds = yt.load(fn_ramses_info)
    # fov = float(ds.length_unit.to('pc').value) * side      # in pc
    unit_l = read_quant_from_ramses_info(fn_ramses_info, 'unit_l')
    boxlen = read_quant_from_ramses_info(fn_ramses_info, 'boxlen')
    fov = boxlen * unit_l * side         # in cm
    half_fov_pc = fov * units.cm / 2

    # if version == "2020":
    if 1:
        ski.setIntAttribute(
            "MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument",
            "numPixelsX", nx)
        ski.setIntAttribute(
            "MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument",
            "numPixelsY", ny)
        ski.setQuantityAttribute(
            "MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument",
            "fieldOfViewX", fov * units.cm, skirtUnit='pc')
        ski.setQuantityAttribute(
            "MonteCarloSimulation/instrumentSystem/InstrumentSystem/instruments/FullInstrument",
            "fieldOfViewY", fov * units.cm, skirtUnit='pc')

    if version == "2020":
        for axis in ['X', 'Y', 'Z']:
            for sign, maxmin in zip([-1, 1], ['min', 'max']):
                ski.setQuantityAttribute("MonteCarloSimulation/mediumSystem/MediumSystem/media/AdaptiveMeshMedium", f"{maxmin}{axis}", sign * half_fov_pc, skirtUnit='pc')
                # ski.setQuantityAttribute("MonteCarloSimulation/mediumSystem/MediumSystem/grid/PolicyTreeSpatialGrid", f"{maxmin}{axis}", sign * half_fov_pc, skirtUnit='pc')
        
    ski.saveTo(skifn)
    logging.info("Set nx = ny = {}".format(ny))
    logging.info("Saved " + skifn)
    return


if __name__ == "__main__":

    # main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])

    # mod_ski(fi=sys.argv[1], fo=sys.argv[2], nph=sys.argv[3], lmax=sys.argv[4],
    #         fn_par=sys.argv[5], fn_nml=sys.argv[6], fn_ramses_info=sys.argv[7])

    print("yes")
