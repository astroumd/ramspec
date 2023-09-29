#!/usr/bin/env python
"""run_amr2cube.py
Create grid data from RAMSES outputs using amr2cube.f90

"""
import os
import subprocess
import f90nml

Thidir = os.path.dirname(os.path.realpath(__file__))
# global EXE
EXE = f"{Thidir}/../tools/amr2cube/amr2cube_mod"


def run_amr2cube_base(jobdir, outs, out_dir, center, width, lma, fields):
    """Run amr2cube.f90 for a given jobdir and output numbers

    Args:
        jobdir (str): the path to the RAMSES job directory
        outs (list of integers):  list of output numbers
        out_dir (str): the path to store the output data
        center (str or list of integers): the center of the sample box. If 'c', then the center is the center of the whole simulation box. Otherwise, it should be a list of three numbers between 0 and 1.
        width (_type_): width of the sample box as a fraction of the simulation box size
        lma (int): the maximum grid levels for grid sampling. The size of the grids would be 1/2^lma of the box size.
        fields (dictionary): a dictionary to store the field names and their corresponding indices
    """

    if isinstance(center, str):
        if center == 'c':
            center = [.5, .5, .5]
        else:
            raise ValueError(f"center={center} is not understood.")
    left = [c - width / 2 for c in center]
    right = [c + width / 2 for c in center]
    params = (f"-xmi {left[0]} -xma {right[0]} -ymi {left[1]} -yma {right[1]} "
              f"-zmi {left[1]} -zma {right[2]}")
    for i in outs:
        print("Running amr2cube for output", i)
        fields_ = ['den', 'xHII']
        indices = [fields[field] for field in fields_]
        for field, typ in zip(fields_, indices):
            denname = f"{out_dir}/out{i:05d}_{field}_l{lma}.dat"
            if os.path.isfile(denname):
                print(f"{denname} exists. Skipping")
                continue
            cmd = (f"{EXE} -inp {jobdir}/output_{i:05d} -out {denname} -typ {typ} -lma {lma} {params}")
            try: 
                ret = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, timeout=300, universal_newlines=True)
            except subprocess.CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
            else:
                print("amr2cube run successfully")

            print(f"{denname} created")


def run_amr2cube_from_nml(ramjobdir, outs, out_dir, nml, fields, lmax=None):
    """Run amr2cube.f90 for a given jobdir and output numbers

    Args:
        ramjobdir (str): the path to the RAMSES job directory
        outs (list of integers): list of output numbers
        out_dir (str): the path to store the output data
        nml (str): the path 
        fields (dictionary): a dictionary to store the field names and their corresponding indices
        lmax (int, optional): the maximum grid levels for grid sampling. The size of the grids would be 1/2^lmax of the box size. Defaults to None, in which case lmax is read from the namelist file.
    """

    params = f90nml.read(nml)["PARAMS"]
    center = [(params["xmax"] + params["xmin"])/2,
                (params["ymax"] + params["ymin"])/2,
                (params["zmax"] + params["zmin"])/2]
    widthx = float(params["xmax"]) - float(params["xmin"])
    widthy = float(params["ymax"]) - float(params["ymin"])
    widthz = float(params["zmax"]) - float(params["zmin"])
    assert(widthx == widthy == widthz, 
           "Non-cubic box is not supported. To fix this, either change it to cubic box or change the code here to define a projection axis.")
    width = widthx

    if lmax is None:
        lmax = int(params["lmax"])

    run_amr2cube_base(ramjobdir, outs, out_dir, center, width, lmax, fields)
    return


# if __name__ == "__main__":

#     def run_amr2cube_of_job(jobid, outs, center='c', width=0.8, lma=9, suffix='', out_dir="../data/data_amr2cube"):

#         sam_dir = "/startrek/chongchong/Sam"
#         os.system(f"mkdir -p {out_dir}/Job{jobid}")
#         out_dir = f"{out_dir}/Job{jobid}"
#         run_amr2cube(f"{sam_dir}/Job{jobid}", outs, out_dir, center=center, width=width, lma=lma, suffix=suffix)
#         return

#     run_amr2cube_of_job('2.2.2', range(15, 49+1))
#     run_amr2cube_of_job('3.2.2', range(14, 44+1, 2))
#     run_amr2cube_of_job('4.2.1', range(14, 48+1, 1))

#     #run_amr2cube('3.2.2', [38], width=0.98, suffix="_w")
#     #run_amr2cube('3.2.2', [26], width=0.98, suffix="_w")
