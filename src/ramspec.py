import os
import numpy as np
import matplotlib.pyplot as plt
import json
from astropy.io import fits
from ..external.ramtools.ramtools import ramses
from .amr2cube import run_amr2cube_from_nml
from .skirt import run_ramski, plot_skirt_spec
from .run_final_spectrum import make_spec_new, combine_ha_with_skirt, make_combined_spec, plot_combined_spec
from .h_alpha import calculate_halpha
from .sed_tools import get_sed, combine_spec


global FIELDS 
FIELDS = {
    'den': 1,
    'xHII': 12,
}


def run_cube_and_ramski(ramjobdir, outs, data_dir, plot_dir, nml, lmax, letdie=True, version="2023"):

    skirt_data_dir = os.path.join(data_dir, "skirt")
    cube_data_dir = os.path.join(data_dir, "amr2cube")
    ha_data_dir = os.path.join(data_dir, "halpha")
    os.makedirs(skirt_data_dir, exist_ok=True)
    os.makedirs(cube_data_dir, exist_ok=True)
    os.makedirs(ha_data_dir, exist_ok=True)
    
    print("\n>>>>>>>>>> Preparing cube data using amr2cube...")
    run_amr2cube_from_nml(ramjobdir, outs, out_dir=cube_data_dir, nml=nml, fields=FIELDS, lmax=lmax)
    print("Success!")

    print("\n>>>>>>>>>> Calculating H-alpha...")
    calculate_halpha(ramjobdir, outs, cube_data_dir, ha_data_dir, lmax, nml=nml, cube_data_index_formater="{:05d}")
    print("Success!")

    print("\n>>>>>>>>>> Running RAMSKI to prepare data for SKIRT...")
    ret = run_ramski(ramjobdir, outs, skirt_data_dir, nml, lmax=lmax, letdie=letdie, version=version)
    if ret == 0:
        print("Success!")
        print(f"Now you can run SKIRT simulations in the output folder: {skirt_data_dir}")
    else:
        print("Failed!")
        return

    print("\n>>>>>>>>>> Optional: plot SKIRT spectrum...")
    plot_skirt_spec(skirt_data_dir, outs, plot_dir)

    feature = "main_mod_total"

    print("\n>>>>>>>>>> Combine SKIRT spectrum with h-alpha...")
    # TODO: split into make_combined_spec and plot_combined_spec
    make_combined_spec(data_dir, outs, feature=feature)
    plot_combined_spec(data_dir, outs, plot_dir, feature=feature)

    return

    print("\n>>>>>>>>>> Post-process SKIRT spectrum...")
    # data_name = "data-v6-cont-die"
    # plot_base_dir = f"../figures/{data_name}"
    # ha_data = "../data/data_amr2cube" # this is v6
    # overwrite = 0
    # feature = "main_CK_corrected_1e6ph_sedsb"
    # boxw = 0.8
    # task = 'all'

    make_spec_new(ramjobdir, skirt_data_dir, plot_base_dir=plot_dir, ha_data=cube_data_dir, fn_nml=nml, outs=outs, fn_fig="halpha", lmax=lmax)

    return
