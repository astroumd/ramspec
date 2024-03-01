from ramspec import ramspec

# make sure to set the following variable for your version of RAMSES
ramspec.FIELDS = {
    'den': 1,
    'xHII': 12,
}


def main():

    # ----- Set the following variables for your application -----
    jobpath = "/Users/cche/Academics/Projects/2017-RAMSES-jobs/Job-full-snaps/Job2.2.2"  # the path to the RAMSES job
    nml = "./ramses2023_test.nml"  # the path to the nml file for SKIRT
    outs = [30] # the output numbers to process
    lmax = 11   # the maximum refinement level for amr2cube (for calculating H-alpha). If you pick half of the simulation box (xmin=0.25, xmax=0.75) in the nml file and set lmax=11 here, then you will get 2^11 * 0.5 = 1024 pixels in each dimension for the amr2cube data and the photograph will have 1024^2 pixels. Note that this is different from the lmax in the nml file, which is for SKIRT radiative transfer calculation and should be set to the highest refinement level in the RAMSES simulation (otherwise SKIRT run will fail, from my experience).
    # ----- End of user-defined variables -----

    data_dir = "./out"   # the base path to store the processed data
    plot_dir = "./out"
    # data_dir = f"{data_base_dir}/test_run_v2023"    # the path to store the processed data
    # plot_dir = f"{data_base_dir}/test_run_v2023"

    ramspec.run_cube_and_ramski(jobpath, outs, data_dir, plot_dir, nml, lmax=lmax)

    return
    

if __name__ == "__main__":

    main()