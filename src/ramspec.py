import os
from .amr2cube import run_amr2cube
from .skirt import run_ramski

global FIELDS 
FIELDS = {
    'den': 1,
    'xHII': 12,
}

def runjob(jobpath, outs, datadir, center, width, lma, nml):
    print(">>>>>>>>>> Preparing cube data using amr2cube...")
    cube_data_path = f"{datadir}/amr2cube"
    os.makedirs(cube_data_path, exist_ok=True)
    run_amr2cube(jobpath, outs, cube_data_path, center, width, lma, fields=FIELDS)

    skirt_data_path = f"{datadir}/skirt"
    run_ramski(jobpath, outs, skirt_data_path, nml, family="CK", letdie=False,
               version="2023")
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

