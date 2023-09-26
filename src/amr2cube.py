#!/usr/bin/env python
"""run_amr2cube.py
Create grid data from RAMSES outputs using amr2cube.f90

"""
import os
import subprocess
import f90nml
from ..external.ramtools.ramtools import ramses

Thidir = os.path.dirname(os.path.realpath(__file__))
# global EXE
EXE = f"{Thidir}/../tools/amr2cube/amr2cube_mod"


def run_amr2cube_base(jobdir, outs, out_dir, center, width, lma, fields):
    """ Run amr2cube.f90 for a given jobdir and output numbers
    """

    # if not os.path.isdir(out_dir):
    #     os.makedirs(out_dir)
        
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
            # process = subprocess.run(cmd, shell=1, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,  universal_newlines=True)
            # print(process.stdout)
            try: 
                output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, timeout=300, universal_newlines=True)
            except subprocess.CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
            else:
                # print("Output: \n{}\n".format(output))
                print("amr2cube run successfully")

            print(f"{denname} created")


def run_amr2cube_from_nml(ramjobdir, outs, out_dir, nml, fields, lmax=None):

    ram = ramses.Ramses(ramjobdir)
    # parse outs
    if isinstance(outs, str):
        if outs == "all":
            outs = ram.get_all_outputs()
        else:
            raise ValueError(f"outs={outs} is not understood.")

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
