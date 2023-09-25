"""run_job_step1.py

Should run on startrek.
1. Create part_CK for all outputs of a job
2. Run RAMSKI to create hydro for all outputs

"""

import os
import sys
import subprocess
import f90nml
from glob import glob
from ..external.ramtools.ramtools import ramses, to_skirt
from .mod_ski import mod_ski

Thidir = os.path.dirname(os.path.realpath(__file__))
# RAMSKI = f"{Thidir}/../external/RAMSKI_test/f90/ramski"
RAMSKI = f"{Thidir}/../external/RAMSKI_test/ramski"
# AMR2CUBE = "/startrek2nb/chongchong/Sam/test_amr2cube/amr2cube_mod"


def run_ramski(ramjobdir, outs, skirt_job_dir, nml, ski=None, family="CK", letdie=True, 
            #    cube_data_dir="out/cube", 
            #    lmax=9, 
               version="2023"):
    """New 2023.4.17: prepare SKIRT simulations for all outputs of a job."""

    os.makedirs(skirt_job_dir, exist_ok=True)

    part = f"part_{family}"
    ram = ramses.Ramses(ramjobdir)
    # parse outs
    if isinstance(outs, str):
        if outs == "all":
            outs = ram.get_all_outputs()
        else:
            raise ValueError(f"outs={outs} is not understood.")

    # nmls = glob(f"{skirt_job_dir}/*.nml")
    # assert len(nmls) == 1, f"More than one nml file in {skirt_job_dir}"
    # fn_nml = nmls[0]
    # ski_name = the_nml["PARAMS"]["name_ski"] + ".ski"
    
    the_nml = f90nml.read(nml)
    params = the_nml["PARAMS"]
    msg = ("{xi}max and {xi}min does not center around 0.5. This is not"
           "consistant with my to_skirt.py file, so I will stop here.")
    iscenter = True
    for xi in ['x', 'y', 'z']:
        iscenter = iscenter and (abs(float(params[xi + "max"]) + float(params[xi + "min"]) - 1) \
                < 1e-10)
    if iscenter:
        center = None   # center is the box center
    else:
        center = [(params["xmax"] + params["xmin"])/2,
                  (params["ymax"] + params["ymin"])/2,
                  (params["zmax"] + params["zmin"])/2]
    width = float(params[xi + "max"]) - float(params[xi + "min"])
    lmax = int(params["lmax"])

    for out in outs:
        print(f"\nRunning RAMSKI for output {out}")
        fn_sink = ram.get_sink_path(out)
        if os.stat(fn_sink).st_size == 0:
            print(fn_sink, 'is empty. Skipped.')
            continue
        
        # make part_CK
        out_dir = f"{skirt_job_dir}/out{out:05d}"
        os.makedirs(out_dir, exist_ok=True)
        fn_part = os.path.join(out_dir, part)
        print(f"Creating {fn_part}")
        to_skirt.to_skirt(jobdir=ramjobdir,
            output=out,
            fn_out=fn_part,
            center=center,
            family=family,
            width_boxlen=width,
            letdie=letdie,
            skip_exist=True,
        )

        # make hydro using RAMSKI
        hydro_fp = os.path.join(out_dir, the_nml["PARAMS"]["name_hdr"])
        inp = ram.get_output_path(out)
        if os.path.isfile(hydro_fp):
            print(f"{out_dir}/hydro exists. Skipped")
        else:
            cmd = f"{RAMSKI} -inp {inp} -nmlpath {nml} -outdir {out_dir}"
            print("running")
            print(cmd)
            try: 
                output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True, timeout=300, universal_newlines=True)
            except subprocess.CalledProcessError as exc:
                print("Status : FAIL", exc.returncode, exc.output)
            else:
                print("RAMSKI runs successfully")
                print(f"Done creating hydro")

        # modify ski file
        if ski is None:
            ski_fn = os.path.join(out_dir, "main.ski")
            ski_mod_fn = ski_fn.replace(".ski", "_mod.ski")
        else:
            ski_fn = ski
            ski_mod_fn = os.path.join(out_dir, "main.ski")

        if os.path.exists(ski_mod_fn):
            print(f"{ski_mod_fn} exists. Skipped.")
        else:
            nph = '1e6'
            fn_info = ram.get_info_path(out)
            mod_ski(ski_fn, ski_mod_fn, nph, lmax, part, nml, fn_info, version=version)

    # # create amr2cube data for halpha calculation
    # # amr2cube_dir = "../data/data_amr2cube/v7"
    # # amr2cube_dir = "1_data/external/tmp/amr2cube_2023"
    # run_amr2cube.EXE = AMR2CUBE
    # run_amr2cube.run_amr2cube_base(ramjobdir, outs, out_dir=cube_data_dir, width=width, lma=lmax)

    return 0


if __name__ == "__main__":

    if len(sys.argv) >= 2:
        run_job(sys.argv[1])
    else:
        print(f"Usage: {sys.argv[0]} jobid")
