import os
import numpy as np
import matplotlib.pyplot as plt
from .amr2cube import run_amr2cube_from_nml
from .skirt import run_ramski, plot_skirt_spec
from .sed_tools import get_sed, combine_spec
from .run_final_spectrum import make_combined_spec, plot_combined_spec
from .h_alpha import calculate_halpha


global FIELDS 
FIELDS = {
    'den': 1,
    'xHII': 12,
}


def run_cube_and_ramski(ramjobdir, outs, data_dir, plot_dir, nml, lmax, letdie=True, version="2023"):
    """Prepare data for SKIRT simulations and calculate H-alpha. When running for a second time after SKIRT simulations are done, will combine SKIRT spectrum with H-alpha and plot the combined spectrum.

    Args:
        ramjobdir (str): the RAMSES job directory (where there should be output_xxxxx folders)
        outs (list): list of integers, the output numbers to process
        data_dir (str): the base path to store the processed data
        plot_dir (str): the path to store the plots
        nml (str): the path to the nml file for SKIRT
        lmax (int): 
        letdie (bool, optional): _description_. Defaults to True.
        version (str, optional): _description_. Defaults to "2023".
    """

    skirt_data_dir = os.path.join(data_dir, "skirt")
    cube_data_dir = os.path.join(data_dir, "amr2cube")
    ha_data_dir = os.path.join(data_dir, "halpha")
    os.makedirs(skirt_data_dir, exist_ok=True)
    os.makedirs(cube_data_dir, exist_ok=True)
    os.makedirs(ha_data_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)
    
    print("\n>>>>>>>>>> Preparing cube data using amr2cube...")
    run_amr2cube_from_nml(ramjobdir, outs, out_dir=cube_data_dir, nml=nml, fields=FIELDS, lmax=lmax)
    print("Success!")

    print("\n>>>>>>>>>> Calculating H-alpha...")
    calculate_halpha(ramjobdir, outs, cube_data_dir, ha_data_dir, lmax, nml=nml)
    print("Success!")

    print("\n>>>>>>>>>> Running RAMSKI to prepare data for SKIRT...")
    ret = run_ramski(ramjobdir, outs, skirt_data_dir, nml, lmax=lmax, letdie=letdie, version=version)
    if ret == 1:
        print("Success!")
        print(f"Now run SKIRT simulations in the output folder: {skirt_data_dir}.")
        print("Then run this script again to continue.")
        return
    elif ret == 2:
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

    # make_spec_new(ramjobdir, skirt_data_dir, plot_base_dir=plot_dir, ha_data=cube_data_dir, fn_nml=nml, outs=outs, fn_fig="halpha", lmax=lmax)

    return
