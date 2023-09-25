import sys, os 
import numpy as np
import matplotlib.pyplot as plt
from academicpython import plottools as pt
from ramtools.ramses import Ramses
# from ramtools.center import get_tff

np.seterr(divide = 'ignore')
LAM_HALPHA = 0.65628  # um
RAMSES = os.path.expanduser("~/Sam")


def broaden_ha(ha, l, width, writeback=False):
    """ Artificially broaden spectrum at ha (without changing the height) """
    lam_l, lam_r = LAM_HALPHA - width / 2, LAM_HALPHA + width / 2
    pick_ha = np.argmax(l >= LAM_HALPHA)
    pick_l = np.argmax(l >= lam_l)
    pick_r = np.argmax(l >= lam_r)
    assert pick_l>0 and pick_r>0, "MyError: failed to find the location of H-alpha in l"
    if pick_r == pick_l:
        pick_r += 1
    if not writeback:
        ret = ha.copy()
        ret[pick_l:pick_r] = ha[pick_ha]
        return ret
    else:
        ha[pick_l:pick_r] = ha[pick_ha]


def plot_sed_ind(jobid, task="combined", letdie=False, broaden=False, shift=False, in_dir=None, out_dir=None):

    # pt.set_plotdir(f"{out_dir}/Job{jobid}")
    pt.set_plotdir(out_dir)
    morestr = "letdie-" if letdie else ""
    # in_dir = f"3_output/data/copied/local_run_v6/Job{jobid}/sed-{morestr}R1000"
    if broaden:
        morestr += "broadenha-"
    if shift:
        morestr += "shift-"
    # figdir = f"data/copied/local_run_v6/Job{jobid}/fig-sed-{morestr}R1000"
    # pt.set_plotdir(figdir)
    # Plot combined spectra
    outs = []
    for i in range(100):
        fn_sed = f"{in_dir}/out{i}.txt"
        if os.path.exists(fn_sed):
            outs.append(i)
    if len(outs) <= 1:
        print(f"Skipped job{jobid}")
        return
    colors = plt.cm.viridis
    plt.figure()
    skip = 2
    count = 0
    for i in outs:
        if i % skip > 0:
            continue
        thecolor = colors((i - outs[0]) / (outs[-1] - outs[0]))
        fn_sed = f"{in_dir}/out{i}.txt"
        l, fll = np.loadtxt(fn_sed, skiprows=2).T
        if broaden:
            broaden_ha(fll, l, width=LAM_HALPHA*(10**0.05-1), writeback=True)
        pick = np.logical_and(l >= 0.01, l <= 1)
        if shift:
            count += 1
            # shift curve to the left
            l *= 10.**(count * -0.02)
        plt.plot(np.log10(l)[pick], np.log10(fll)[pick], color=thecolor)
    plt.xlabel(r"$\log \lambda$ [$\mu$m]")
    plt.ylabel(r"$\log (F_\lambda \lambda)$ [erg/s/cm$^2$]")
    # plt.ylim([1e-12, 1e-8])
    # plt.xlim([10**-2.1, 10**0.1])
    plt.ylim([-12, -6])
    plt.xlim([-2.1, 0.1])
    # plt.grid(True)
    # pt.save(f"sed-combined-job{jobid}", isprint=1)
    pt.save(f"sed-combined-{morestr}job{jobid}", isprint=1)
    if task == 'combined':
        return

    # Plot individual figures
    for i in range(100):
        fn_sed = f"{in_dir}/out{i}.txt"
        fo = f"fig-out{i}"
        if not os.path.exists(fn_sed):
            continue
        l, fll = np.loadtxt(fn_sed, skiprows=2).T
        pick = np.logical_and(l >= 0.01, l <= 1)
        plt.figure()
        plt.loglog(l[pick], fll[pick])
        plt.xlabel(r"$\lambda$ [$\mu$m]")
        plt.ylabel(r"$F_\lambda \lambda$ [erg/s/cm$^2$]")
        plt.ylim([1e-11, 1e-7])
        plt.xlim([10**-2.1, 10**0.1])
        plt.grid(True)
        figfn = f"fig-out{i}"
        # if sty == 'talk':
        #     figfn = 'talk-' + figfn
        fn = pt.save(figfn)
        plt.close()
