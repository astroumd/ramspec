#!/startrek/chongchong/anaconda3/envs/yt/bin/python -u

import os
from glob import glob
from src.process import run_job
from src.run_amr2cube import run_amr2cube

Out_fp_fmt = "/startrek2nb/chongchong/Sam/Job{jobid}/output_00001/info_00001.txt"

def do_a_job(jobid):

    # v6
    skirt_job_dir="../data/yorp07/run_v6"
    part = "part_CK_corrected"
    name = "main_CK_corrected"
    skip = 2
    width = 0.8
    amr2cube_dir = "../data/data_amr2cube"

    # v7
    # skirt_job_dir = "../data/yorp07/run_v7"
    # part = "part"
    # name = "main"
    # skip = 4
    # width = 0.98
    # amr2cube_dir = "../data/data_amr2cube/v7"

    # {{{ make hydro and part
    run_job(jobid, outs=None, skip=skip, skirt_job_dir=skirt_job_dir,
            part=part, letdie=True)
    run_job(jobid, outs=None, skip=skip, skirt_job_dir=skirt_job_dir,
            part=part+"_notdie", letdie=False)
    # }}}

    # {{{ mode ski
    nml = f"main-job{jobid}.nml"
    #ski_fn = nml.replace(".nml", ".ski")
    ski_fn = glob(f"{skirt_job_dir}/Job{jobid}/out*")[0] + "/main.ski"
    ski_mod_fn = f"{skirt_job_dir}/Job{jobid}/{name}.ski"
    ski_mod_fn = ski_mod_fn.replace(".ski", "_1e6ph.ski")
    ski_nd_mod_fn = f"{skirt_job_dir}/Job{jobid}/{name}_nd.ski"
    ski_nd_mod_fn = ski_nd_mod_fn.replace(".ski", "_1e6ph.ski")
    #skirt_job_dir_injob = os.path.join(skirt_job_dir, "Job" + jobid)
    for ski, partn in zip([ski_mod_fn, ski_nd_mod_fn], [part, part+"_notdie"]):
        if os.path.exists(ski):
            print(f"{ski} exists. Skipped modifying ski file")
        else:
            # old
            #mod_ski(...)
            cmd = "{cwd}/mod_ski.py {fi} {fo} {nph} {lmax} {fn_par} {skirt_job_dir}/{fn_nml} {fn_ramses_info}".format(
                skirt_job_dir=skirt_job_dir,
                cwd='.',
                fi=ski_fn,
                fo=ski,
                nph='1e6',
                lmax=9,
                fn_par=partn,
                fn_nml=nml,
                fn_ramses_info=Out_fp_fmt.format(jobid=jobid),
            )
            print(cmd)
            os.system(cmd)
    # }}}

    # {{{
    #sed_out_fmt = (f"{skirt_job_dir}/Job{jobid}/out{{out:02d}}/"
    #               f"{name}_1e6ph_sedsb_sed.dat")
    for i in range(1, 100):
        if not os.path.exists(f"{skirt_job_dir}/Job{jobid}/out{i:02d}"):
            continue
        if not i % skip == 0:
            continue
        run_amr2cube(jobid, [i], out_dir=amr2cube_dir, width=width)
    # }}}

    print("Now, run run_all_CK.sh in skirt job directory.")

    return

def main():

    # jobs = ['4.2.1', '2.2.2', '3.2.2']
    jobs = ['2.2.2']
    for jobid in jobs:
        do_a_job(jobid)

main()
