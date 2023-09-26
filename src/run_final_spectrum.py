#!/startrek/chongchong/anaconda3/envs/yt/bin/python -u

""" Make final spectrum, combining continuum with H-alpha """

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
# import astropy.units as u, astropy.constants as c
import json
import f90nml

from academicpython import plottools as pt
from . import h_alpha
from ..external.ramtools.ramtools import utilities, ramses

matplotlib.rcParams['figure.dpi'] = 300
plt.style.use(['science', 'no-latex'])

C_hu_halpha = 3.026e-12  # in erg, = 1.89 eV
C_arcsec2_to_sr = 2.35044e-11
lam_halpha = 0.65628  # um

def scaled(x):
    return x / x.max()

def mass_to_lum(m):
    """ Convert mass (Msun) to luminosity (Lsun). Source: wikipedia
    mass-luminosity relation
    """

    if m < 0.43:
        return 0.23 * m**2.3
    if m < 2:
        return m**4.0
    if m < 55:
        return 1.4 * m**3.5
    else:
        return 32000 * m

def set_y_decades(decades, ax=None):
    if ax is None:
        ax = plt.gca()
    ymin = ax.get_ylim()[1] / 10**decades
    ax.set_ylim(bottom=ymin)

def main_plot_compare_ratio():
    """ Compare line ratio (halpha / continuum) """
    # plot line ratio v.s. density
    _, ax = pt.scaled_figure()
    jobs = ['4.2.1', '2.2.2', '3.2.2']
    dens = [2e2, 2e3, 2e4]
    plot_base_dir = "../figures/run_v7"
    for den, jobid in zip(dens, jobs):
        t, ratio = make_spec(plot_base_dir, jobid, is_ret_ratio=True)
        ax.plot(t, ratio, label="n = {:.1e} cm$^{{-3}}$".format(den))
    ax.set(xlabel=r"$t/t_{\rm ff}$", ylabel=r"$F_{H\alpha}/F_{\rm cont.}$")
    ax.legend(loc=0)
    pt.set_plotdir(plot_base_dir)
    pt.save("Joball-line-ratio", fromfile=__file__)

def main2():
    plot_base_dir = "../results/fig/run_v6"
    f, axs = plt.subplots(3, 3, figsize=[12,6], sharex='col', )
    #context = ['science']
    jobids = ['4.2.1', '2.2.2', '3.2.2']
    #with plt.style.context(context):
    if 1:
        for i in range(3):
            make_spec(plot_base_dir, jobids[i], is_overwrite_json=0, axs=axs[:, i])
        for axrow in axs:
            for axcol in axrow:
                axcol.get_legend().remove()
        for axrow in axs[:, -1]:
            axrow.legend(bbox_to_anchor=(1.01, 0.5), loc='center left',
                         fontsize='medium')
        for axi in axs[-1, :]:
            axi.set_ylim([0, 1.1])
        pt.clean_sharey_label(axs)
        plt.subplots_adjust(hspace=0.02)
        pt.set_plotdir(plot_base_dir)
        pt.save("three_jobs", f)



def make_spec_new(ram_job_dir, data_job_dir, plot_base_dir, ha_data, fn_nml, outs, is_overwrite_json=0, axs=None, is_ret_ratio=False, feature=None, fn_fig="ha-{jobid}", data_name="", R=None, is_times_lambda=False, lmax=9, task='all'):

    r = ramses.Ramses(ram_job_dir)
    halpha_data_base_dir = ha_data
    the_nml = f90nml.read(fn_nml)
    params = the_nml["PARAMS"]
    box_fraction = float(params["xmax"]) - float(params["xmin"])
    boxw = box_fraction
    fn_json = f"{plot_base_dir}/data.json"
    side = box_fraction * r.unit_l  # cm
    #ha_width = 0.6 / 8  # micron, JWST, source at z=8
    ha_width = lam_halpha / 100  # l/dl = 100
    if feature is None:
        feature = "{}_mod_total".format(params["name_ski"])

    if (not os.path.exists(fn_json)) or is_overwrite_json:
        h_alpha.do_full_spec(
            ram_job_dir,
            data_job_dir,
            plot_base_dir,
            halpha_data_base_dir,
            fn_json,
            outs,
            box_fraction,
            feature=feature,
            ha_width=ha_width,
            lmax=lmax,
            cubeformater="{:05d}",
        )
    else:
        print("json file exists. Not going to redo.")
    return


def make_spec_new_cont():
    """ Continuation of make_spec_new """

    jobname = os.path.basename(ram_job_dir)

    with open(fn_json, 'r') as fin:
        data = json.load(fin)
    if task == 'sed':
        return

    # ef = ef_qesc.EfPlot(jobid, jobdir=f"{HOME}/Sam")
    ef = ef_qesc.EfPlot(jobpath=ram_job_dir)
    ef.params['mass_shift'] = 0.3
    ef.params['is_rescale_Q'] = True
    ef.params['is_always_rescale_Q'] = True
    Q = ef.get_Q_tot_time_list(4, 'HI', params=ef.params)
    QEsc = ef.get_Qesc_tot_time_list(4, 'HI', params=ef.params)
    fesc = QEsc / Q
    out_fesc_end = len(fesc)

    outs = [i for i in data['outs'] if i <= out_fesc_end]
    data_cut = len(outs)
    halpha_str = np.array(data['halpha_strength'])[:data_cut]
    # unit: flux per arcsec2
    halpha_no_dust_str = np.array(data['halpha_no_dust_strength'])[:data_cut]
    halpha_height = halpha_str / data['halpha_width']
    halpha_no_dust_height = halpha_no_dust_str / data['halpha_width']
    str_cont = []
    str_cont_otherl = []
    lam_other = 1.0
    for out in outs:
        pick = np.argmax(np.array(data['dust_wavelength'][str(out)]) >= lam_halpha)
        str_cont.append(data['dust_spec'][str(out)][pick])
        pick_otherl = np.argmax(np.array(data['dust_wavelength'][str(out)])
                                >= lam_other)
        str_cont_otherl.append(data['dust_spec'][str(out)][pick_otherl])
    str_cont = np.array(str_cont)  # (len(outs), ) array
    str_cont_otherl = np.array(str_cont_otherl)
    ratio = halpha_height / str_cont  # (len(outs), ) array
    # compare halpha strength with (1 - fesc) * m_*
    # ax.plot(outs, (list_ha_str / width) / halpha_dust_str, label="H-alpha / continuum")
    #outs_fesc = np.arange(out_fesc_end) + 1 # np.array(outs)[np.array(outs) <= out_fesc_end]
    #outs_fesc = np.array(outs)[np.array(outs) <= out_fesc_end]
    outs_fesc = outs
    Q_outs = Q[np.array(outs_fesc) - 1]
    f_trap_outs = np.array([1 - fesc[out-1] for out in outs_fesc])
    times = [r.get_time(out) - r.tRelax for out in outs]
    mtl = np.vectorize(mass_to_lum)
    mtt = np.vectorize(utilities.mass_to_temp)
    mtr = np.vectorize(utilities.mass_to_radius)
    times_fesc = []
    mstars = []
    mstars_die = []
    lum_stars_die = []
    Tsum_stars_die = []
    rayleis = []
    Q_notdie_scaled = []
    Q_inside_not_die = []
    Q_notdie = []
    for out in outs_fesc:
        times_fesc.append(r.get_time(out) - r.tRelax)
        m = r.get_sink_masses(out)
        alive = r.is_sink_alive(out)
        mstars.append(m.sum())
        mstars_die.append(np.sum(m[alive]))
        lum_stars_die.append(np.sum(mtl(m[alive] * 0.4)))
        Tsum_stars_die.append(np.sum(mtt(m[alive] * 0.4)))
        raylei = np.sum(mtr(m[alive] * 0.4)**2 * mtt(m[alive] * 0.4))
        rayleis.append(raylei)
        Q_notdie_scaled.append(np.sum(ef.get_Q(out, let_star_die=False)))
        # 0.450 * hv * Q
        # eV = 1.60218e-12                            # erg
        # # hv_ha = 13.6 * (1/2**2 - 1/3**2) * eV       # erg
        # Q_inside_not_die.append(ef.get_Q(out, let_star_die=0).sum())
        Q_inside_not_die.append(
            ef.get_Q_in_box(out, let_star_die=0, boxw=boxw).sum() * unit_lum)
        Q_notdie.append(np.sum(ef.get_Q(out, let_star_die=False)) * unit_lum)

    # ndot = Q_outs * f_trap_outs  # unit_lum s^-1
    ndot = Q_inside_not_die * f_trap_outs  # unit_lum s^-1
    # hu_halpha = 5.447e-12  # in erg, = 3.4eV
    hu_halpha = C_hu_halpha                # erg
    ha_predict = 0.450 * hu_halpha * ndot  # erg s^-1
    ha_str_predict = ha_predict * 1. / (4 * np.pi * side**2) # erg s^-1 cm^-2 sr^-1
    ha_str_predict *= C_arcsec2_to_sr # (n_outs,) erg s-1 cm-2 arcsec-2

    times_fesc = np.array(times_fesc)
    mstars = np.array(mstars)
    mstars_die = np.array(mstars_die)
    lum_stars_die = np.array(lum_stars_die)
    Tsum_stars_die = np.array(Tsum_stars_die)
    rayleis = np.array(rayleis)
    Q_notdie_scaled = np.array(Q_notdie_scaled)
    Q_inside_not_die = np.array(Q_inside_not_die)
    Q_notdie = np.array(Q_notdie)
    #outs_fesc = np.array(outs)[np.array(outs) <= out_fesc_end]
    f_trap_times_mstar = f_trap_outs * mstars
    f_trap_times_mstar_die = f_trap_outs * mstars_die
    f_trap_times_Q_die = f_trap_outs * Q_outs
    f_trap_times_Q = f_trap_outs * Q_notdie_scaled

    if is_ret_ratio:
        assert len(times_fesc) == len(ratio), \
                "{} {}".format(len(times_fesc), len(ratio))
        # return t/tff, halpha_height / continuum
        return times_fesc / ef.tff, ratio

    # jobid_simple = jobid.replace('.', '')
    jobid_simple = jobname
    fn_comb = f"{plot_base_dir}/{jobname}/combined_spec_v4-{jobid_simple}.pdf"
    if not os.path.exists(fn_comb):
        colors = plt.cm.viridis
        _, ax = plt.subplots()
        # Spectrums
        #for i in outs[::2]:
        for i in outs:
            color = colors((i - outs[0]) / (outs[-1] - outs[0]))
            ax.plot(data['dust_wavelength'][str(i)],
                    np.log10(data['combined_spec'][str(i)]), color=color)
        ax.set(xscale='log', xlabel=r"$\lambda$ [$\mu$m]", yscale='linear',
               ylabel=r'$\log F_{\lambda}$ [erg/(s cm2 arcsec2 $\mu$m)]',
               xlim=[1e-2, 2e0],
               #ylim=[6e-18, 3e-11],
        )
        #_ = pt.set_y_decades(6, ax, is_ret_ymax=1)
        pt.set_y_height(8.3, ax, is_ret_ymax=0)
        pt.setMyLogFormat(ax, 'x')
        #pt.set_plotdir(f"{plot_base_dir}/Job{jobid}")
        #pt.save(f"combined_spec_v4.pdf")
        plt.savefig(fn_comb)
        plt.close()

    is_save = True
    if axs is None:
        _, axs = plt.subplots(4, 1, figsize=[5,8])
        context = ['science']
    else:
        context = []
        is_save = False

    def plot_rayleigh(ax):
        # Normalized quantity: continuum, Rayleigh-Jeans
        #ax.plot(times_fesc, Q_outs / np.max(Q_outs), '-', label=r"$Q_{\rm cl}$")
        # ax.plot(times, str_cont / np.max(str_cont), '-',
        #         label=r"Continuum at H$\alpha$")
        ax.plot(times_fesc, scaled(Q_notdie), '-', label=r"$L$(LyC)")
        ax.plot(times, scaled(str_cont), '-', label=r"Continuum at H$\alpha$")
        # ax.plot(times, scaled(str_cont_otherl), '-',
        #         label=r"Continuum at {:.2f} um".format(lam_other))
        # ax.plot(times_fesc, lum_stars_die / np.max(lum_stars_die),
        #         '-', label=r"$L_{\rm cl}$")
        ax.plot(times_fesc, rayleis / np.max(rayleis), '-',
                label=r"Rayleigh-Jeans spectrum at H$\alpha$")
        # ax.plot(times_fesc, mstars_die / np.max(mstars_die), '-',
        #         label=r"$m_{\rm \star, die}$")
        # Q(Halpha)
        # ax.set_ylabel("Normalized quantity")
        ax.set_yscale("log")
        ax.legend(loc=0, fontsize='small')
        # pt.set_y_decades(3, ax=ax)

    def plot_ha_cont(ax):
        # Ha vs continuum
        ax.plot(times, np.log10(halpha_height), 'C0-', label=r"H$\alpha$",)
        # ax.plot(times, halpha_no_dust_height, 'C0--',
        #         label=r"H-$\alpha$, no dust")
        ax.plot(times, np.log10(str_cont), 'C1-', label="Continuum")
        # ax.set_ylabel(r"$I$ [erg/(s cm2 arcsec2 $\mu$m)]", fontsize='small')
        yl = r"$I$ (erg s-1 cm-2 arcsec-2 $\mu$-1)"
        yl = "log " + yl
        ax.set_ylabel(yl, fontsize='small')
        # ax.set_yscale('log')
        # pt.set_y_decades(3, ax=ax)
        ax.legend(loc=0, fontsize='small')
        #pt.save("h_alpha_strength.pdf")
        #plt.show()

        ##ax = axs[1]
        ##ax.plot(times, np.array(str_cont), 'C1-', label=r"continuum")
        ###ax.set_ylabel(r"$F_{\rm H\alpha} / F_{\rm continuum}$ at H$\alpha$")
        ###ax.legend(bbox_to_anchor=(0.5, 1.15), )
        ##ax.legend(loc='upper left')
        ###pt.save("halpha_to_dust.pdf")

        #ax = axs[2]
        #ax.plot(times, halpha_height / np.array(str_cont),
        #        label=r"$\frac{{\rm H}\alpha}{\rm continuum}$")
        #       #label=r"H-$\alpha$ to continuum ratio")
        ##ax.set_ylabel(r"$F_{\rm H\alpha} / F_{\rm continuum}$ at H$\alpha$")
        ##ax.legend(bbox_to_anchor=(0.5, 1.15), )
        #ax.set_ylim(bottom=0)
        #ax.legend(loc='upper right')
        ##pt.save("halpha_to_dust.pdf")

        # integrated spectrum
        # f, ax = plt.subplots()
        # ax.plot(times, [inte(dust_spec[out], dust_wavelength[out]) for out in outs])
        # plt.show()
        return

    def plot_norm_ha(ax):
        # normalized quantity: Q(1 - fesc), Ha
        ax.plot(times, halpha_str / halpha_str.max(), 'C0-', label=r"H$_{\rm alpha}$")
        ax.plot(times, halpha_no_dust_str / halpha_no_dust_str.max(), 'C0--',
                label=r"H$_{\rm alpha}$, no dust")
        # ax.plot(times_fesc, f_trap_times_mstar / f_trap_times_mstar.max(), 'C1-',
        #         label=r"$m_{\star} (1 - f_{\rm esc})$")
        #ax.plot(times_fesc, f_trap_times_mstar_die / f_trap_times_mstar_die.max(), 'C1--',
        #        label=r"$m_{\rm \star, alive} (1 - f_{\rm esc})$")
        #ratio = halpha_str.max() / f_trap_times_mstar_die.max()
        ax.plot(times_fesc, f_trap_times_Q / f_trap_times_Q.max(), 'C2-',
                label=r"$Q_{\rm \star} (1 - f_{\rm esc})$")
        ratio = halpha_str.max() / f_trap_times_Q.max()
        # ax.text(.98, .02, "ratio = {:.2e}".format(ratio), va='bottom', ha='right',
        #         transform=ax.transAxes)
        ax.set_ylim(0, 1.3)
        ax.set_ylabel("Normalized quantity", fontsize='small')
        #ax.legend(bbox_to_anchor=(0.5, 1.25), ncol=3, loc='upper center', fontsize='small')
        ax.legend(ncol=2, loc='upper right', fontsize='x-small')
        ax.set(xlabel='t (Myr)', )

    def plot_norm_ha_simple(ax):    # 2022-03-14
        ax.plot(times, halpha_str / halpha_str.max(), 'C0-', label=r"H$_{\rm alpha}$")
        # ax.plot(times, halpha_no_dust_str / halpha_no_dust_str.max(), 'C0--',
        #         label=r"H$_{\rm alpha}$, no dust")
        ax.plot(times_fesc, f_trap_times_Q / f_trap_times_Q.max(), 'C2-',
                label=r"$Q_{\rm \star} (1 - f_{\rm esc})$")
        ratio = halpha_str.max() / f_trap_times_Q.max()
        ax.set_ylim(0, 1.3)
        ax.set_ylabel("Normalized quantity", fontsize='small')
        ax.legend(ncol=2, loc='upper right', fontsize='x-small')
        ax.set(xlabel='t (Myr)', )

    def plot_fesc_and_Q(ax):
        ax.plot(times_fesc, f_trap_outs, label="1 - fesc")
        # ax.plot(times_fesc, Q_notdie_scaled / np.max(Q_notdie_scaled), label=r"$Q_\star$")
        # Qscale = np.max(Q_inside_not_die)
        Qscale = 1e50
        ax.plot(times_fesc, Q_inside_not_die / Qscale,
                label=r"$Q_{{\rm *, in box}}$ / {:.1e}".format(Qscale))
        ax.plot(times_fesc, Q_notdie / Qscale,
                label=r"$Q_*$ / {:.1e}".format(Qscale))
        ax.set(yscale='log', ylim=[1e-3, 20])
        ax.legend(loc=0)

    def plot_flux(ax, is_take_log):
        d_nd = halpha_no_dust_str if not is_take_log else np.log10(halpha_no_dust_str)
        d_d = halpha_str if not is_take_log else np.log10(halpha_str)
        d_m = ha_str_predict if not is_take_log else np.log10(ha_str_predict)
        ax.plot(times, d_nd, 'C0--', label="data (no dust)")
        ax.plot(times, d_nd, 'C0.', )
        # ax.plot(times, d_d, 'C1--', label="data")
        # ax.plot(times, d_d, 'C1.', )
        ax.plot(times_fesc, d_m, 'C2-', label="model")
        ax.legend(loc=0)
        yl = r"$F_{\rm H\alpha}$ (erg s-1 cm-2 arcsec-2)"
        if is_take_log:
            yl = "log " + yl
        ax.set_ylabel(yl, fontsize='small')

    def plot_ratio(ax):
        # plot ratio of model to data
        ax.plot(times_fesc, ha_str_predict / halpha_no_dust_str[:len(ha_str_predict )],
                'C0-', )
        ax.set(xlabel="t (Myr)", ylim=[0, 4],
               ylabel="model / data",
              )
        #ax.legend(loc='upper right')

    def plot_ha_to_cont(ax):
        ax.plot(times, halpha_height / str_cont, )
        ax.set(ylim=[1e0, 1e2], yscale='log')
        ax.set_ylabel(r"$F(H\alpha) / F({\rm cont})$")

    def plot_complex():
        # _, axs = plt.subplots(4, 1, figsize=[5,8])
        _, axs = pt.scaled_figure(8, 1)

        plot_norm_ha(axs[0])
        plot_flux(axs[1], True)
        plot_flux(axs[2], False)
        plot_ratio(axs[3])
        plot_fesc_and_Q(axs[4])
        plot_ha_cont(axs[5])
        plot_rayleigh(axs[6])
        plot_ha_to_cont(axs[7])

        pt.clean_sharex(axs)
        plt.subplots_adjust(hspace=0.02)
        pt.set_plotdir(plot_base_dir)
        #pt.save(f"h-alpha-Job{jobid}")
        pt.save(fn_fig.format(jobid=jobid), fromfile=__file__)

    def plot_single():

        for is_take_log in [True, False]:
            _, ax = pt.scaled_figure()
            datay = halpha_no_dust_str if not is_take_log else np.log10(halpha_no_dust_str)
            datay1 = halpha_str if not is_take_log else np.log10(halpha_str)
            datay2 = ha_str_predict if not is_take_log else np.log10(ha_str_predict)
            ax.plot(times, datay, 'C0--', label="data (no dust)")
            ax.plot(times, datay, 'C0.', )
            ax.plot(times, datay1, 'C1--', label="data")
            ax.plot(times, datay1, 'C1.', )
            ax.plot(times_fesc, datay2, 'C2-', label="model")
            ax.legend(loc=0)
            # ax.set_ylabel("Halpha strength (erg s-1 cm-2 arcsec-1)", fontsize='x-small')
            yl = r"$F_{\rm H\alpha}$ (erg s-1 cm-2 arcsec-2)"
            if is_take_log:
                yl = "log " + yl
                fn = f"ha-fesc-Job{jobid}.pdf"
            else:
                #ax.set_yscale('log')
                fn = f"ha-fesc-lin-Job{jobid}.pdf"
            ax.set_ylabel(yl, fontsize='small')
            pt.set_plotdir(plot_base_dir)
            pt.save(fn)

    def plot_normed():
        _, ax = pt.scaled_figure()
        plot_norm_ha(ax)
        pt.set_plotdir(plot_base_dir)
        pt.save(fn_fig.format(jobid=jobid) + "-normed")

        _, ax = pt.scaled_figure()
        plot_norm_ha_simple(ax)
        pt.set_plotdir(plot_base_dir)
        pt.save(fn_fig.format(jobid=jobid) + "-normed-simple")

    #plot_complex()
    # plot_fig2()
    # plot_single()
    plot_normed()
    
    return


def make_spec(data_base_dir, plot_base_dir, ha_data, jobid=None, is_overwrite_json=0, axs=None, is_ret_ratio=False, feature=None, fn_fig="ha-{jobid}", boxw=1.0, data_name="", R=None, is_times_lambda=False, task='all', ram_job_dir=None):

    # data_job_dir = data_base_dir
    if ram_job_dir is None:
        assert jobid is not None
        ram_job_dir = f"{HOME}/Sam/Job{jobid}"

    data_job_dir = f"{data_base_dir}/Job{jobid}"

    make_spec_new(ram_job_dir, data_base_dir, plot_base_dir, ha_data, is_overwrite_json=is_overwrite_json, axs=axs, is_ret_ratio=is_ret_ratio, feature=feature, fn_fig=fn_fig, boxw=boxw, data_name=data_name, R=R, is_times_lambda=is_times_lambda, task=task)
    return


def make_spec_new2(ram_job_dir, data_job_dir, plot_base_dir, ha_data, is_overwrite_json=0, axs=None, is_ret_ratio=False, feature=None, fn_fig="ha-{jobid}", boxw=1.0, data_name="", R=None, is_times_lambda=False, task='all'):


    # fn_nml = f"{data_base_dir}/main-job{jobid}.nml"
    the_nml = f90nml.read(fn_nml)
    params = the_nml["PARAMS"]
    box_fraction = float(params["xmax"]) - float(params["xmin"])

    fn_json = f"{plot_base_dir}/Job{jobid}/data.json"
    # side = box_fraction * r.boxlen * r.unit_l  # cm. WRONG!!!!
    side = box_fraction * r.unit_l  # cm

    is_doit = 1
    if os.path.exists(fn_json):
        if type(is_overwrite_json) is str:
            if is_overwrite_json == 'bk':
                pt.backup(fn_json)
        elif not is_overwrite_json:
            print("json file exists. Not going to redo")
            is_doit = 0
    #ha_width = 0.6 / 8  # micron, JWST, source at z=8
    ha_width = lam_halpha / 100  # l/dl = 100
    if is_doit:
        h_alpha.do_full_spec(
            jobid,
            data_job_dir,
            plot_base_dir,
            halpha_data_base_dir,
            fn_json,
            outs=None,
            feature=feature,
            box_fraction=box_fraction,
            ha_width=ha_width,
        )
    with open(fn_json, 'r') as fin:
        data = json.load(fin)
    if task == 'sed':
        return

    ef = ef_qesc.EfPlot(jobid, jobdir=f"{HOME}/Sam")
    ef.params['mass_shift'] = 0.3
    ef.params['is_rescale_Q'] = True
    ef.params['is_always_rescale_Q'] = True
    Q = ef.get_Q_tot_time_list(4, 'HI', params=ef.params)
    QEsc = ef.get_Qesc_tot_time_list(4, 'HI', params=ef.params)
    fesc = QEsc / Q
    out_fesc_end = len(fesc)

    outs = [i for i in data['outs'] if i <= out_fesc_end]
    data_cut = len(outs)
    halpha_str = np.array(data['halpha_strength'])[:data_cut]
    # unit: flux per arcsec2
    halpha_no_dust_str = np.array(data['halpha_no_dust_strength'])[:data_cut]
    halpha_height = halpha_str / data['halpha_width']
    halpha_no_dust_height = halpha_no_dust_str / data['halpha_width']
    str_cont = []
    str_cont_otherl = []
    lam_other = 1.0
    for out in outs:
        pick = np.argmax(np.array(data['dust_wavelength'][str(out)]) >= lam_halpha)
        str_cont.append(data['dust_spec'][str(out)][pick])
        pick_otherl = np.argmax(np.array(data['dust_wavelength'][str(out)])
                                >= lam_other)
        str_cont_otherl.append(data['dust_spec'][str(out)][pick_otherl])
    str_cont = np.array(str_cont)  # (len(outs), ) array
    str_cont_otherl = np.array(str_cont_otherl)
    ratio = halpha_height / str_cont  # (len(outs), ) array
    # compare halpha strength with (1 - fesc) * m_*
    # ax.plot(outs, (list_ha_str / width) / halpha_dust_str, label="H-alpha / continuum")
    #outs_fesc = np.arange(out_fesc_end) + 1 # np.array(outs)[np.array(outs) <= out_fesc_end]
    #outs_fesc = np.array(outs)[np.array(outs) <= out_fesc_end]
    outs_fesc = outs
    Q_outs = Q[np.array(outs_fesc) - 1]
    f_trap_outs = np.array([1 - fesc[out-1] for out in outs_fesc])
    times = [r.get_time(out) - r.tRelax for out in outs]
    mtl = np.vectorize(mass_to_lum)
    mtt = np.vectorize(utilities.mass_to_temp)
    mtr = np.vectorize(utilities.mass_to_radius)
    times_fesc = []
    mstars = []
    mstars_die = []
    lum_stars_die = []
    Tsum_stars_die = []
    rayleis = []
    Q_notdie_scaled = []
    Q_inside_not_die = []
    Q_notdie = []
    for out in outs_fesc:
        times_fesc.append(r.get_time(out) - r.tRelax)
        m = r.get_sink_masses(out)
        alive = r.is_sink_alive(out)
        mstars.append(m.sum())
        mstars_die.append(np.sum(m[alive]))
        lum_stars_die.append(np.sum(mtl(m[alive] * 0.4)))
        Tsum_stars_die.append(np.sum(mtt(m[alive] * 0.4)))
        raylei = np.sum(mtr(m[alive] * 0.4)**2 * mtt(m[alive] * 0.4))
        rayleis.append(raylei)
        Q_notdie_scaled.append(np.sum(ef.get_Q(out, let_star_die=False)))
        # 0.450 * hv * Q
        # eV = 1.60218e-12                            # erg
        # # hv_ha = 13.6 * (1/2**2 - 1/3**2) * eV       # erg
        # Q_inside_not_die.append(ef.get_Q(out, let_star_die=0).sum())
        Q_inside_not_die.append(
            ef.get_Q_in_box(out, let_star_die=0, boxw=boxw).sum() * unit_lum)
        Q_notdie.append(np.sum(ef.get_Q(out, let_star_die=False)) * unit_lum)

    # ndot = Q_outs * f_trap_outs  # unit_lum s^-1
    ndot = Q_inside_not_die * f_trap_outs  # unit_lum s^-1
    # hu_halpha = 5.447e-12  # in erg, = 3.4eV
    hu_halpha = C_hu_halpha                # erg
    ha_predict = 0.450 * hu_halpha * ndot  # erg s^-1
    ha_str_predict = ha_predict * 1. / (4 * np.pi * side**2) # erg s^-1 cm^-2 sr^-1
    ha_str_predict *= C_arcsec2_to_sr # (n_outs,) erg s-1 cm-2 arcsec-2

    times_fesc = np.array(times_fesc)
    mstars = np.array(mstars)
    mstars_die = np.array(mstars_die)
    lum_stars_die = np.array(lum_stars_die)
    Tsum_stars_die = np.array(Tsum_stars_die)
    rayleis = np.array(rayleis)
    Q_notdie_scaled = np.array(Q_notdie_scaled)
    Q_inside_not_die = np.array(Q_inside_not_die)
    Q_notdie = np.array(Q_notdie)
    #outs_fesc = np.array(outs)[np.array(outs) <= out_fesc_end]
    f_trap_times_mstar = f_trap_outs * mstars
    f_trap_times_mstar_die = f_trap_outs * mstars_die
    f_trap_times_Q_die = f_trap_outs * Q_outs
    f_trap_times_Q = f_trap_outs * Q_notdie_scaled

    if is_ret_ratio:
        assert len(times_fesc) == len(ratio), \
                "{} {}".format(len(times_fesc), len(ratio))
        # return t/tff, halpha_height / continuum
        return times_fesc / ef.tff, ratio

    jobid_simple = jobid.replace('.', '')
    fn_comb = f"{plot_base_dir}/Job{jobid}/combined_spec_v4-{jobid_simple}.pdf"
    if not os.path.exists(fn_comb):
        colors = plt.cm.viridis
        _, ax = plt.subplots()
        # Spectrums
        #for i in outs[::2]:
        for i in outs:
            color = colors((i - outs[0]) / (outs[-1] - outs[0]))
            ax.plot(data['dust_wavelength'][str(i)],
                    np.log10(data['combined_spec'][str(i)]), color=color)
        ax.set(xscale='log', xlabel=r"$\lambda$ [$\mu$m]", yscale='linear',
               ylabel=r'$\log F_{\lambda}$ [erg/(s cm2 arcsec2 $\mu$m)]',
               xlim=[1e-2, 2e0],
               #ylim=[6e-18, 3e-11],
        )
        #_ = pt.set_y_decades(6, ax, is_ret_ymax=1)
        pt.set_y_height(8.3, ax, is_ret_ymax=0)
        pt.setMyLogFormat(ax, 'x')
        #pt.set_plotdir(f"{plot_base_dir}/Job{jobid}")
        #pt.save(f"combined_spec_v4.pdf")
        plt.savefig(fn_comb)
        plt.close()

    is_save = True
    if axs is None:
        _, axs = plt.subplots(4, 1, figsize=[5,8])
        context = ['science']
    else:
        context = []
        is_save = False

    def plot_rayleigh(ax):
        # Normalized quantity: continuum, Rayleigh-Jeans
        #ax.plot(times_fesc, Q_outs / np.max(Q_outs), '-', label=r"$Q_{\rm cl}$")
        # ax.plot(times, str_cont / np.max(str_cont), '-',
        #         label=r"Continuum at H$\alpha$")
        ax.plot(times_fesc, scaled(Q_notdie), '-', label=r"$L$(LyC)")
        ax.plot(times, scaled(str_cont), '-', label=r"Continuum at H$\alpha$")
        # ax.plot(times, scaled(str_cont_otherl), '-',
        #         label=r"Continuum at {:.2f} um".format(lam_other))
        # ax.plot(times_fesc, lum_stars_die / np.max(lum_stars_die),
        #         '-', label=r"$L_{\rm cl}$")
        ax.plot(times_fesc, rayleis / np.max(rayleis), '-',
                label=r"Rayleigh-Jeans spectrum at H$\alpha$")
        # ax.plot(times_fesc, mstars_die / np.max(mstars_die), '-',
        #         label=r"$m_{\rm \star, die}$")
        # Q(Halpha)
        # ax.set_ylabel("Normalized quantity")
        ax.set_yscale("log")
        ax.legend(loc=0, fontsize='small')
        # pt.set_y_decades(3, ax=ax)

    def plot_ha_cont(ax):
        # Ha vs continuum
        ax.plot(times, np.log10(halpha_height), 'C0-', label=r"H$\alpha$",)
        # ax.plot(times, halpha_no_dust_height, 'C0--',
        #         label=r"H-$\alpha$, no dust")
        ax.plot(times, np.log10(str_cont), 'C1-', label="Continuum")
        # ax.set_ylabel(r"$I$ [erg/(s cm2 arcsec2 $\mu$m)]", fontsize='small')
        yl = r"$I$ (erg s-1 cm-2 arcsec-2 $\mu$-1)"
        yl = "log " + yl
        ax.set_ylabel(yl, fontsize='small')
        # ax.set_yscale('log')
        # pt.set_y_decades(3, ax=ax)
        ax.legend(loc=0, fontsize='small')
        #pt.save("h_alpha_strength.pdf")
        #plt.show()

    ##ax = axs[1]
    ##ax.plot(times, np.array(str_cont), 'C1-', label=r"continuum")
    ###ax.set_ylabel(r"$F_{\rm H\alpha} / F_{\rm continuum}$ at H$\alpha$")
    ###ax.legend(bbox_to_anchor=(0.5, 1.15), )
    ##ax.legend(loc='upper left')
    ###pt.save("halpha_to_dust.pdf")

    #ax = axs[2]
    #ax.plot(times, halpha_height / np.array(str_cont),
    #        label=r"$\frac{{\rm H}\alpha}{\rm continuum}$")
    #       #label=r"H-$\alpha$ to continuum ratio")
    ##ax.set_ylabel(r"$F_{\rm H\alpha} / F_{\rm continuum}$ at H$\alpha$")
    ##ax.legend(bbox_to_anchor=(0.5, 1.15), )
    #ax.set_ylim(bottom=0)
    #ax.legend(loc='upper right')
    ##pt.save("halpha_to_dust.pdf")

    # integrated spectrum
    # f, ax = plt.subplots()
    # ax.plot(times, [inte(dust_spec[out], dust_wavelength[out]) for out in outs])
    # plt.show()

    def plot_norm_ha(ax):
        # normalized quantity: Q(1 - fesc), Ha
        ax.plot(times, halpha_str / halpha_str.max(), 'C0-', label=r"H$_{\rm alpha}$")
        ax.plot(times, halpha_no_dust_str / halpha_no_dust_str.max(), 'C0--',
                label=r"H$_{\rm alpha}$, no dust")
        # ax.plot(times_fesc, f_trap_times_mstar / f_trap_times_mstar.max(), 'C1-',
        #         label=r"$m_{\star} (1 - f_{\rm esc})$")
        #ax.plot(times_fesc, f_trap_times_mstar_die / f_trap_times_mstar_die.max(), 'C1--',
        #        label=r"$m_{\rm \star, alive} (1 - f_{\rm esc})$")
        #ratio = halpha_str.max() / f_trap_times_mstar_die.max()
        ax.plot(times_fesc, f_trap_times_Q / f_trap_times_Q.max(), 'C2-',
                label=r"$Q_{\rm \star} (1 - f_{\rm esc})$")
        ratio = halpha_str.max() / f_trap_times_Q.max()
        # ax.text(.98, .02, "ratio = {:.2e}".format(ratio), va='bottom', ha='right',
        #         transform=ax.transAxes)
        ax.set_ylim(0, 1.3)
        ax.set_ylabel("Normalized quantity", fontsize='small')
        #ax.legend(bbox_to_anchor=(0.5, 1.25), ncol=3, loc='upper center', fontsize='small')
        ax.legend(ncol=2, loc='upper right', fontsize='x-small')
        ax.set(xlabel='t (Myr)', )

    def plot_norm_ha_simple(ax):    # 2022-03-14
        ax.plot(times, halpha_str / halpha_str.max(), 'C0-', label=r"H$_{\rm alpha}$")
        # ax.plot(times, halpha_no_dust_str / halpha_no_dust_str.max(), 'C0--',
        #         label=r"H$_{\rm alpha}$, no dust")
        ax.plot(times_fesc, f_trap_times_Q / f_trap_times_Q.max(), 'C2-',
                label=r"$Q_{\rm \star} (1 - f_{\rm esc})$")
        ratio = halpha_str.max() / f_trap_times_Q.max()
        ax.set_ylim(0, 1.3)
        ax.set_ylabel("Normalized quantity", fontsize='small')
        ax.legend(ncol=2, loc='upper right', fontsize='x-small')
        ax.set(xlabel='t (Myr)', )

    # _, axs = plt.subplots(4, 1, figsize=[5, 12])

    def plot_fesc_and_Q(ax):
        ax.plot(times_fesc, f_trap_outs, label="1 - fesc")
        # ax.plot(times_fesc, Q_notdie_scaled / np.max(Q_notdie_scaled), label=r"$Q_\star$")
        # Qscale = np.max(Q_inside_not_die)
        Qscale = 1e50
        ax.plot(times_fesc, Q_inside_not_die / Qscale,
                label=r"$Q_{{\rm *, in box}}$ / {:.1e}".format(Qscale))
        ax.plot(times_fesc, Q_notdie / Qscale,
                label=r"$Q_*$ / {:.1e}".format(Qscale))
        ax.set(yscale='log', ylim=[1e-3, 20])
        ax.legend(loc=0)

    # # ax = axs[1]
    # # ax.plot(times_fesc, Q_notdie_scaled, label="Q")
    # # ax.plot(times_fesc, Q_inside_not_die, label="Q inside")
    # # ax.legend(loc=0)

    def plot_flux(ax, is_take_log):
        d_nd = halpha_no_dust_str if not is_take_log else np.log10(halpha_no_dust_str)
        d_d = halpha_str if not is_take_log else np.log10(halpha_str)
        d_m = ha_str_predict if not is_take_log else np.log10(ha_str_predict)
        ax.plot(times, d_nd, 'C0--', label="data (no dust)")
        ax.plot(times, d_nd, 'C0.', )
        # ax.plot(times, d_d, 'C1--', label="data")
        # ax.plot(times, d_d, 'C1.', )
        ax.plot(times_fesc, d_m, 'C2-', label="model")
        ax.legend(loc=0)
        yl = r"$F_{\rm H\alpha}$ (erg s-1 cm-2 arcsec-2)"
        if is_take_log:
            yl = "log " + yl
        ax.set_ylabel(yl, fontsize='small')

    def plot_ratio(ax):
        # plot ratio of model to data
        ax.plot(times_fesc, ha_str_predict / halpha_no_dust_str[:len(ha_str_predict )],
                'C0-', )
        ax.set(xlabel="t (Myr)", ylim=[0, 4],
               ylabel="model / data",
              )
        #ax.legend(loc='upper right')

    def plot_ha_to_cont(ax):
        ax.plot(times, halpha_height / str_cont, )
        ax.set(ylim=[1e0, 1e2], yscale='log')
        ax.set_ylabel(r"$F(H\alpha) / F({\rm cont})$")

    def plot_complex():
        # _, axs = plt.subplots(4, 1, figsize=[5,8])
        _, axs = pt.scaled_figure(8, 1)

        plot_norm_ha(axs[0])
        plot_flux(axs[1], True)
        plot_flux(axs[2], False)
        plot_ratio(axs[3])
        plot_fesc_and_Q(axs[4])
        plot_ha_cont(axs[5])
        plot_rayleigh(axs[6])
        plot_ha_to_cont(axs[7])

        pt.clean_sharex(axs)
        plt.subplots_adjust(hspace=0.02)
        pt.set_plotdir(plot_base_dir)
        #pt.save(f"h-alpha-Job{jobid}")
        pt.save(fn_fig.format(jobid=jobid), fromfile=__file__)

    def plot_single():

        for is_take_log in [True, False]:
            _, ax = pt.scaled_figure()
            datay = halpha_no_dust_str if not is_take_log else np.log10(halpha_no_dust_str)
            datay1 = halpha_str if not is_take_log else np.log10(halpha_str)
            datay2 = ha_str_predict if not is_take_log else np.log10(ha_str_predict)
            ax.plot(times, datay, 'C0--', label="data (no dust)")
            ax.plot(times, datay, 'C0.', )
            ax.plot(times, datay1, 'C1--', label="data")
            ax.plot(times, datay1, 'C1.', )
            ax.plot(times_fesc, datay2, 'C2-', label="model")
            ax.legend(loc=0)
            # ax.set_ylabel("Halpha strength (erg s-1 cm-2 arcsec-1)", fontsize='x-small')
            yl = r"$F_{\rm H\alpha}$ (erg s-1 cm-2 arcsec-2)"
            if is_take_log:
                yl = "log " + yl
                fn = f"ha-fesc-Job{jobid}.pdf"
            else:
                #ax.set_yscale('log')
                fn = f"ha-fesc-lin-Job{jobid}.pdf"
            ax.set_ylabel(yl, fontsize='small')
            pt.set_plotdir(plot_base_dir)
            pt.save(fn)

    def plot_normed():
        _, ax = pt.scaled_figure()
        plot_norm_ha(ax)
        pt.set_plotdir(plot_base_dir)
        pt.save(fn_fig.format(jobid=jobid) + "-normed")

        _, ax = pt.scaled_figure()
        plot_norm_ha_simple(ax)
        pt.set_plotdir(plot_base_dir)
        pt.save(fn_fig.format(jobid=jobid) + "-normed-simple")


    #plot_complex()
    # plot_fig2()
    # plot_single()
    plot_normed()

def combine_specs():
    # combine H-alpha with continuum

    fn_nml = f"{data_base_dir}/main-job{jobid}.nml"
    the_nml = f90nml.read(fn_nml)
    params = the_nml["PARAMS"]
    box_fraction = float(params["xmax"]) - float(params["xmin"])
    r = ramses.Ramses(jobid=jobid, jobdir=f"{HOME}/Sam")
    side = box_fraction * r.unit_l  # cm

    lam_halpha = 0.65628   # um
    #ha_width = 0.6 / 8  # micron, JWST, source at z=8
    ha_width = lam_halpha / 100  # l/dl = 100
    is_doit = 1

def main():

    # data_name = "data-v6-cont-notdie"
    # data_base_dir = f"../data/yorp07/run_v6"
    # plot_base_dir = "../figures/run_v6"
    # ha_data = "../data/data_amr2cube" # this is v6
    # overwrite = 0
    # feature = "main_CK_corrected_nd_1e6ph_sedsb"
    # fn_fig = "ha-{jobid}"
    # boxw = 0.8

    #data_name = "data-v6-cont-die"
    #data_base_dir = f"../data/yorp07/run_v6"
    #plot_base_dir = f"../figures/run_v6"
    #ha_data = "../data/data_amr2cube" # this is v6
    #overwrite = 0
    ##feature = "main_CK_corrected_nd_1e6ph_sedsb"
    #feature = "main_CK_corrected_1e6ph_sedsb"
    #fn_fig = "ha-{jobid}"
    #boxw = 0.8

    data_name = "data-v6-cont-die"
    data_base_dir = f"../data/yorp07/run_v6"
    plot_base_dir = f"../figures/{data_name}"
    ha_data = "../data/data_amr2cube" # this is v6
    overwrite = 0
    feature = "main_CK_corrected_1e6ph_sedsb" # {data_base_dir}/Job{jobid}/{feature}.ski
    boxw = 0.8
    task = 'all'

    # # 2022-03-10, rerun in 2022. Add the following
    task = 'sed'
    overwrite = 1

    #data_base_dir = f"../data/yorp07/run_v7"
    #plot_base_dir = "../figures/run_v7"
    #ha_data = "../data/data_amr2cube/v7"
    #feature = None

    # jobids = ['4.2.1', '2.2.2', '3.2.2']
    # jobids = ['4.2.1']
    jobids = ['2.2.2']
    for jobid in jobids:
        make_spec(data_base_dir, plot_base_dir, ha_data, jobid,
                  is_overwrite_json=overwrite, feature=feature,
                  boxw=boxw, data_name=data_name, task=task)
    return


if __name__ == "__main__":
    main()
    #main_plot_compare_ratio()
    #main2()
    #tmp1()
