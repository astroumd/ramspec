"""run_job_step1.py

Should run on startrek.
1. Create part_CK for all outputs of a job
2. Run RAMSKI to create hydro for all outputs

"""

import os
import subprocess
import f90nml
from glob import glob
from ..external.ramtools.ramtools import ramses, to_skirt
from .mod_ski import mod_ski


Thidir = os.path.dirname(os.path.realpath(__file__))
RAMSKI = f"{Thidir}/../external/RAMSKI_test/f90/ramski"


def run_ramski(ramjobdir, outs, skirt_job_dir, nml, lmax=None, 
               letdie=True, version="2023"):
    """New 2023.4.17: prepare SKIRT simulations for all outputs of a job."""

    os.makedirs(skirt_job_dir, exist_ok=True)

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

    SEDtype = f90nml.read(nml)["SOURCE"]["SEDtype"]
    part = f"part_{SEDtype}"

    ski_base_name = params["name_ski"]

    if lmax is None:
        lmax = params["lmax"]

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
            family=SEDtype,
            width_boxlen=width,
            letdie=letdie,
            skip_exist=True,
        )

        # make hydro using RAMSKI
        hydro_fp = os.path.join(out_dir, params["name_hdr"])
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
                return 1
            else:
                print("RAMSKI runs successfully")
                print(f"Done creating hydro")

        # modify ski file
        ski_fn = os.path.join(out_dir, ski_base_name + ".ski")
        ski_mod_fn = ski_fn.replace(".ski", "_mod.ski")
        if os.path.exists(ski_mod_fn):
            print(f"{ski_mod_fn} exists. Skipped.")
        else:
            fn_info = ram.get_info_path(out)
            mod_ski(ski_fn, ski_mod_fn, lmax, part, nml, fn_info, version=version)

    return 0


# if __name__ == "__main__":

#     if len(sys.argv) >= 2:
#         run_job(sys.argv[1])
#     else:
#         print(f"Usage: {sys.argv[0]} jobid")
