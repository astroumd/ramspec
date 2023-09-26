import os
import f90nml
from .amr2cube import run_amr2cube_from_nml
from .skirt import run_ramski
from ..external.ramtools.ramtools import ramses
from ..external.ramtools.ramtools.to_skirt import to_skirt
from .run_final_spectrum import make_spec_new


global FIELDS 
FIELDS = {
    'den': 1,
    'xHII': 12,
}


def run_cube_and_ramski(ramjobdir, outs, skirt_data_dir, cube_data_dir, nml, lmax, letdie=True, version="2023"):
    
    print(">>>>>>>>>> Preparing cube data using amr2cube...")
    os.makedirs(cube_data_dir, exist_ok=True)
    run_amr2cube_from_nml(ramjobdir, outs, out_dir=cube_data_dir, nml=nml, fields=FIELDS, lmax=lmax)
    print("Success!")
    print()

    print(">>>>>>>>>> Running RAMSKI...")
    ret = run_ramski(ramjobdir, outs, skirt_data_dir, nml, lmax=lmax, letdie=letdie, version=version)
    if ret == 0:
        print("Success!")
        print()
        print(f"Now you can run SKIRT simulations in the output folder: {skirt_data_dir}")
    else:
        print("Failed!")
        print()

    print(">>>>>>>>>> Calculate H-alpha spectrum...")
    # data_name = "data-v6-cont-die"
    # plot_base_dir = f"../figures/{data_name}"
    # ha_data = "../data/data_amr2cube" # this is v6
    # overwrite = 0
    # feature = "main_CK_corrected_1e6ph_sedsb"
    # boxw = 0.8
    # task = 'all'

    # # 2022-03-10, rerun in 2022. Add the following
    task = 'sed'
    overwrite = 1

    jobids = ['2.2.2']
    for jobid in jobids:
        make_spec_new(ramjobdir, 
                      skirt_data_dir, 
                      plot_base_dir='.',
                      ha_data=cube_data_dir, 
                      fn_nml=nml,
                      outs=outs,
                      fn_fig="halpha",
                      lmax=lmax,
                      )

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

