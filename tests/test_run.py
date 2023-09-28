import os
from ramspec import ramspec

# make sure to set the following variable for your version of RAMSES
ramspec.FIELDS = {
    'den': 1,
    'xHII': 12,
}

# this is the relative/absolute path to store the processed data
DATA_PATH = "data"

def main2020():
    # set the path to the RAMSES job here
    # jobpath = "/Users/cche/Academics/Projects/2017-RAMSES-jobs/Job-full-snaps/Job2.0.v2.full"
    jobpath = "/Users/cche/Academics/Projects/2017-RAMSES-jobs/Job-full-snaps/Job2.2.2"
    datadir = f"{DATA_PATH}/test_run3"
    nml = "./ramses2020.nml"
    outs = [40]
    ramspec.runjob(jobpath, outs, datadir, center='c', width=0.8, lma=9, nml=nml)
    print("Success!")


def main():

    # ----- Set the following variables for your application -----
    data_base_dir = "../data"   # the base path to store the processed data
    jobpath = "/Users/cche/Academics/Projects/2017-RAMSES-jobs/Job-full-snaps/Job2.2.2"  # the path to the RAMSES job
    nml = "../nml/ramses2023_test.nml"  # the path to the nml file
    outs = [30] # the output numbers to process
    lmax = 11   # the maximum refinement level for amr2cube (for calculating H-alpha). Note that this is different from the lmax in the nml file. 
    data_dir = f"{data_base_dir}/test_run_v2023"    # the path to store the processed data
    plot_dir = f"{data_base_dir}/test_run_v2023"

    ramspec.run_cube_and_ramski(jobpath, outs, data_dir, plot_dir, nml, lmax=lmax, letdie=False, version="2023")

    return
    

if __name__ == "__main__":

    main()