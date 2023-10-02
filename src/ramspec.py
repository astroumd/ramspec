import os
import f90nml
from .amr2cube import run_amr2cube_base
from .skirt import run_ramski, plot_skirt_spec
from .sed_tools import get_sed, combine_spec
# from .run_final_spectrum import make_combined_spec, plot_combined_spec, make_mock_photo
from .h_alpha import read_dat, read_dat_v2, calculate_halpha
from .spectrum import *
from .plot_image import plot_halpha, plot_combined_spec


global FIELDS 
FIELDS = {
    'den': 1,
    'xHII': 12,
}


def run_cube_and_ramski(ramjobdir, outs, data_dir, plot_dir, nml, lmax, letdie=False):
    """Prepare data for SKIRT simulations and calculate H-alpha. When running for a second time after SKIRT simulations are done, will combine SKIRT spectrum with H-alpha and plot the combined spectrum.

    Args:
        ramjobdir (str): the RAMSES job directory (where there should be output_xxxxx folders)
        outs (list): list of integers, the output numbers to process
        data_dir (str): the base path to store the processed data
        plot_dir (str): the path to store the plots
        nml (str): the path to the nml file for SKIRT
        lmax (int): 
        letdie (bool, optional): whether or not to allow stars 'die' based on their age. Defaults to True.
    """

    skirt_data_dir = os.path.join(data_dir, "skirt")
    cube_data_dir = os.path.join(data_dir, "amr2cube")
    ha_data_dir = os.path.join(data_dir, "halpha")
    os.makedirs(skirt_data_dir, exist_ok=True)
    os.makedirs(cube_data_dir, exist_ok=True)
    os.makedirs(ha_data_dir, exist_ok=True)
    os.makedirs(plot_dir, exist_ok=True)

    params = f90nml.read(nml)["PARAMS"]
    left = [params["xmin"], params["ymin"], params["zmin"]]
    right = [params["xmax"], params["ymax"], params["zmax"]]
    # center = [(params["xmax"] + params["xmin"])/2,
    #             (params["ymax"] + params["ymin"])/2,
    #             (params["zmax"] + params["zmin"])/2]
    widthx = float(params["xmax"]) - float(params["xmin"])
    widthy = float(params["ymax"]) - float(params["ymin"])
    widthz = float(params["zmax"]) - float(params["zmin"])
    assert widthx == widthy == widthz, "Non-cubic box is not supported"
    width = widthx
    version = params["version"]
    
    print("\n>>>>>>>>>> Preparing cube data using amr2cube...")
    # run_amr2cube_from_nml(ramjobdir, outs, out_dir=cube_data_dir, nml=nml, fields=FIELDS, lmax=lmax)
    run_amr2cube_base(ramjobdir, outs, out_dir=cube_data_dir, left=left, right=right, lmax=lmax, fields=FIELDS) 
    print("Success!")

    print("\n>>>>>>>>>> Calculating H-alpha...")
    calculate_halpha(ramjobdir, outs, cube_data_dir, ha_data_dir, lmax, nml=nml)
    print("Success!")

    print("\n>>>>>>>>>> Plotting H-alpha image...")
    plot_halpha(data_dir, outs, plot_dir, box_fraction=width)

    print("\n>>>>>>>>>> Running RAMSKI to prepare data for SKIRT...")
    ret = run_ramski(ramjobdir, outs, skirt_data_dir, nml, lmax=lmax, letdie=letdie, version=version)
    if ret == 1:
        print("Success!")
        print(f"Now run SKIRT simulations in the following folder: {skirt_data_dir}/out#####. To run a SKIRT simulation, do `skirt main_mod.ski`. You may change the number of photons in the SKIRT simulation by editing the file main_mod.ski file by hand before running SKIRT.")
        print("Then run this script again to continue.")
        return
    elif ret == 2:
        print("Failed!")
        return

    print("\n>>>>>>>>>> Optional: plot SKIRT spectrum...")
    plot_skirt_spec(skirt_data_dir, outs, plot_dir)     # empty for now

    feature = "main_mod_total" 
    # feature = "main_mod_1e7ph_total"

    print("\n>>>>>>>>>> Combine SKIRT spectrum with h-alpha...")
    # TODO: split into make_combined_spec and plot_combined_spec
    make_combined_spec(data_dir, outs, feature=feature)

    print("\n>>>>>>>>>> Plotting combined spectrum...")
    plot_combined_spec(data_dir, outs, plot_dir, feature=feature)

    print("\n>>>>>>>>>> Apply filters to combined data cube...")
    apply_filters(data_dir, outs, feature)

    return
