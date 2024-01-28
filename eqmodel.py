from math import log, log10, log2
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from scipy.optimize import fsolve
import warnings
import mplcursors
from SwitchHelperFunctions import *


class DimericSwitch():
    def __init__(self, kd_unphosphorylated, kd_phosphorylated):
        self.kd_unphosphorylated = kd_unphosphorylated
        self.kd_phosphorylated = kd_phosphorylated

    def __str__(self):
        return str(f'Dimeric switch with a Kd of {self.kd_unphosphorylated} M when unphosphorylated and of {self.kd_phosphorylated} M when fully phosphorylated.')

    def get_key_values(self):
        optimise_concentration(self.kd_unphosphorylated,
                               self.kd_phosphorylated, printresults=True)

    def get_key_values_at_dimer_baseline(self, b=0.05):
        optimise_dimer_baseline(self.kd_unphosphorylated,
                                self.kd_phosphorylated, b, printresults=True)

    def get_key_values_at_monomer_baseline(self, b=0.05):
        optimise_monomer_baseline(
            self.kd_unphosphorylated, self.kd_phosphorylated, b, printresults=True)

    def get_behaviour(self, *args, **kwargs):
        switch_behaviour(self.kd_unphosphorylated,
                         self.kd_phosphorylated, *args, **kwargs)

    def get_kd(self):
        return (self.kd_unphosphorylated, self.kd_phosphorylated)


def dimer_concentration(A0, B0, Kd):
    Sum = A0 + B0 + Kd
    AB = (Sum - np.sqrt(Sum ** 2 - 4 * A0 * B0))/2
    A_free = A0 - AB
    B_free = B0 - AB
    results = {
        'Dimer': AB,
        'A_free': A_free,
        'B_free': B_free,
        'Dimer fraction': AB / min(A0, B0)
    }
    return results


def fret_model(donor, acceptor, Kd, lowerlimit=39, upperlimit=155):
    Sum = donor + acceptor + Kd
    AB = (Sum - np.sqrt(Sum ** 2 - 4 * donor * acceptor))/2
    donor_free = donor - AB
    acceptor_free = acceptor - AB
    results = {
        'Dimer': AB,
        'Free donor': donor_free,
        'Free acceptor': acceptor_free,
        'Bound donor fraction': AB / donor,
        'Bound acceptor fraction': AB / acceptor,
        'Ratiometric FRET': lowerlimit + ((AB / donor) * (upperlimit - lowerlimit))
    }
    return results


def on_off_change(kd_np, kd_p, c, printresults=False):
    warnings.filterwarnings('ignore', message='divide by zero')
    warnings.filterwarnings('ignore', message='invalid value encountered')
    dimers_np = dimer_concentration(c, c, kd_np)['Dimer fraction']
    dimers_p = dimer_concentration(c, c, kd_p)['Dimer fraction']
    difference = dimers_p - dimers_np
    if (-1e-6 < difference < 1e-6) or dimers_np < 0.00001 or dimers_p < 0.00001:
        fold_change = kd_np / kd_p
    else:
        fold_change = dimers_p / dimers_np
    ratio = fold_change if fold_change >= 1 else 1/fold_change
    results = {
        'Fraction difference': difference,
        'Fraction_np': dimers_np,
        'Fraction_p': dimers_p,
        'Fold-change': fold_change,
    }
    if printresults == True:
        c = "{:.3e}".format(c)
        dotSeparation(
            f"Dimer fraction at {c} M (unphosphorylated):", "{:.3f}".format(dimers_np))
        dotSeparation(
            f"Dimer fraction at {c} M (fully phosphorylated):", "{:.3f}".format(dimers_p))
        dotSeparation(f"Difference in dimer fraction:",
                      "{:.3f}".format(difference))
        if fold_change >= 1:
            dotSeparation(
                f"Fold-change in dimer fraction:",
                "{:.3f}".format(fold_change),
                "("+"{:.3f}".format(ratio)+" X more dimers)"
            )
        else:
            dotSeparation(
                f"Fold-change in dimer fraction:",
                "{:.3f}".format(fold_change),
                "("+"{:.3f}".format(ratio)+" X less dimers)"
            )
    return results


def optimise_concentration(kd1, kd2, printresults=False):
    def get_diff(x):
        return on_off_change(kd1, kd2, x)['Fraction difference']
    # Find the closest order of magnitude
    xlist = np.logspace(-14, 0, 15)
    fit = x_for_maxabs_y(xlist, get_diff)
    # Scan around the closest order of magnitude
    xmin, xmax = log10(fit / 10), log10(fit * 10)
    xlist = np.logspace(xmin, xmax, 1000)
    fit = x_for_maxabs_y(xlist, get_diff)
    # Finer scan around previous best fit
    xmin, xmax = log10(fit / 3), log10(fit * 3)
    xlist = np.logspace(xmin, xmax, 10000)
    optimal_concentration = x_for_maxabs_y(xlist, get_diff)
    # Report results
    results = on_off_change(kd1, kd2, optimal_concentration)
    results.update({'Optimal concentration': optimal_concentration})
    if printresults == True:
        dotSeparation(
            f"Greatest absolute change in dimer fraction at:",
            "{:.3e}".format(optimal_concentration),
            "M"
        )
        on_off_change(kd1, kd2, optimal_concentration, printresults=True)
    return results


def optimise_dimer_baseline(kd1, kd2, baseline=0.05, printresults=False):
    def get_dimer_fraction(x):
        a = max(kd1, kd2)
        b = min(kd1, kd2)
        y = on_off_change(a, b, x)['Fraction_np']
        y = y - baseline
        return y
    # Find the closest order of magnitude
    xlist = np.logspace(-14, 0, 15)
    fit = x_for_minabs_y(xlist, get_dimer_fraction)
    # Scan around the closest order of magnitude
    xmin, xmax = log10(fit / 10), log10(fit * 10)
    xlist = np.logspace(xmin, xmax, 1000)
    fit = x_for_minabs_y(xlist, get_dimer_fraction)
    # Finer scan around previous best fit
    xmin, xmax = log10(fit / 3), log10(fit * 3)
    xlist = np.logspace(xmin, xmax, 10000)
    conc_for_baseline = x_for_minabs_y(xlist, get_dimer_fraction)
    # Report results
    results = on_off_change(kd1, kd2, conc_for_baseline)
    results.update({'Concentration for baseline': conc_for_baseline})
    if printresults == True:
        bl = "{:.2f}".format(baseline)
        dotSeparation(
            f"The concentration for a {bl} dimer baseline is:",
            "{:.3e}".format(conc_for_baseline),
            "M"
        )
        on_off_change(kd1, kd2, conc_for_baseline, printresults=True)
    return results


def optimise_monomer_baseline(kd1, kd2, baseline=0.05, printresults=False):
    def get_monomer_fraction(x):
        a = max(kd1, kd2)
        b = min(kd1, kd2)
        y = on_off_change(a, b, x)['Fraction_p']
        y = y - (1 - baseline)
        return y
    # Find the closest order of magnitude
    xlist = np.logspace(-14, 0, 15)
    fit = x_for_minabs_y(xlist, get_monomer_fraction)
    # Scan around the closest order of magnitude
    xmin, xmax = log10(fit / 10), log10(fit * 10)
    xlist = np.logspace(xmin, xmax, 1000)
    fit = x_for_minabs_y(xlist, get_monomer_fraction)
    # Finer scan around previous best fit
    xmin, xmax = log10(fit / 3), log10(fit * 3)
    xlist = np.logspace(xmin, xmax, 10000)
    conc_for_baseline = x_for_minabs_y(xlist, get_monomer_fraction)
    # Report results
    results = on_off_change(kd1, kd2, conc_for_baseline)
    results.update({'Concentration for baseline': conc_for_baseline})
    if printresults == True:
        bl = "{:.2f}".format(baseline)
        dotSeparation(
            f"The concentration for a {bl} monomer baseline is:",
            "{:.3e}".format(conc_for_baseline),
            "M"
        )
        on_off_change(kd1, kd2, conc_for_baseline, printresults=True)
    return results


def behaviour_data(kd1, kd2, xlist):
    ylist1 = []
    ylist2 = []
    ylist3 = []
    ylist4 = []
    for x in xlist:
        y1 = on_off_change(kd1, kd2, x)['Fraction_np']
        y2 = on_off_change(kd1, kd2, x)['Fraction_p']
        y3 = on_off_change(kd1, kd2, x)['Fraction difference']
        y4 = on_off_change(kd1, kd2, x)['Fold-change']
        ylist1.append(y1)
        ylist2.append(y2)
        ylist3.append(y3)
        ylist4.append(y4)
    return (xlist, ylist1, ylist2, ylist3, ylist4)


def switch_behaviour(kd_np, kd_p, *args, type=None, annotations=True):
    if not type in [None, "optimal", "dimerbaseline", "monomerbaseline"]:
        raise ValueError(
            "Graph type must be 'optimal' or 'dimerbaseline', 'monomerbaseline'.")

    if type == None:
        list1 = args
        args = []
        for i in list1:
            if isinstance(i, float):
                args.append(i)
        args.sort()

    if type == "optimal":
        opt = optimise_concentration(kd_np, kd_p)
        c = opt['Optimal concentration']
        args = [c]

    if type == 'dimerbaseline':
        list1 = args
        args = []
        for i in list1:
            b = optimise_dimer_baseline(kd_np, kd_p, i)
            b_conc = b['Concentration for baseline']
            args.append(b_conc)
        args.sort()

    if type == 'monomerbaseline':
        list1 = args
        args = []
        for i in list1:
            b = optimise_monomer_baseline(kd_np, kd_p, i)
            b_conc = b['Concentration for baseline']
            args.append(b_conc)
        args.sort()

    xlist = np.logspace(-14, -1, 1000)
    xlist, ylist1, ylist2, ylist3, ylist4 = behaviour_data(kd_np, kd_p, xlist)

    fig, ax = plt.subplots(1, 2, figsize=(10, 5.4))

    ax[0].set_axisbelow(True)
    ax[0].grid(which='major', axis='y')
    ax[0].plot(xlist, ylist1, lw=2,
               label=f'$K_d$ = {kd_np} M (Unphosphorylated)')
    ax[0].plot(xlist, ylist2, lw=2, label=f'$K_d$ = {kd_p} M (Phosphorylated)')
    ax[0].fill_between(xlist, ylist1, ylist2, color='#EEEEEE', alpha=0.8)
    ax[0].set_xscale('log')
    ax[0].set_xlim(1e-14, 1e-1)
    ax[0].set_xticks([1e-12, 1e-9, 1e-6, 1e-3])
    ax[0].set_ylim(0, 1)
    ax[0].set_title(f'Dimerisation propensity', size=12)
    ax[0].set_xlabel(f'Switch concentration (M)', size=12)
    ax[0].set_ylabel(f'Dimer fraction', size=12)
    ax[0].legend(loc="upper left", bbox_to_anchor=(-0.025, -0.25),
                 ncol=1, edgecolor='none', fontsize=10, labelspacing=0.75)

    ax[1].set_xscale('log')
    ax[1].set_xlim(1e-14, 1e-1)
    ax[1].set_xticks([1e-12, 1e-9, 1e-6, 1e-3])
    ax[1].set_yscale('log')
    ax[1].set_ylim(1e-3, 1e3)
    ax[1].plot(xlist, ylist4, linestyle='--', lw=1, c='green',
               label='Dimer fraction fold-change after phosphorylation (r)')
    ax[1].set_title(f'Effect of phosphorylation on dimerisation', size=12)
    ax[1].set_xlabel(f'Switch concentration (M)', size=12)
    ax[1].set_ylabel(f'Fold-change', size=12)
    ax[1].legend(loc="upper left", bbox_to_anchor=(-0.07, -0.25),
                 ncol=1, edgecolor='none', fontsize=10, labelspacing=0.75)

    ax_delta = ax[1].twinx()
    ax[1].set_zorder(1)
    ax[1].patch.set_alpha(0)
    ax_delta.plot([0, 1], [0, 0], c='black', ls='-', lw=1)
    ax_delta.set_ylim(-1, 1)
    ax_delta.plot(xlist, ylist3, lw=2, c='gray',
                  label='Additional dimer fraction after phosphorylation ($\Delta$)')
    ax_delta.fill_between(xlist, ylist3, 0, color='#EEEEEE', alpha=0.8)
    ax_delta.set_ylabel(f'Difference', size=12)
    ax_delta.legend(loc="upper left", bbox_to_anchor=(-0.07, -0.35),
                    ncol=1, edgecolor='none', fontsize=10, labelspacing=0.75)

    n = 1
    for arg in args:
        r = on_off_change(kd_np, kd_p, arg)
        fraction_np = r['Fraction_np']
        fraction_p = r['Fraction_p']
        ax[0].vlines(arg, -1, max(fraction_np, fraction_p),
                     color='gray', linestyle=':')
        globals()[f'a{n}'] = ax[0].plot(
            [arg, arg],
            [fraction_np, fraction_p],
            marker='o',
            ms=3.5, linestyle='',
            color="black"
        )
        n += 1
        ax[0].plot([arg, arg], [fraction_np, fraction_p], color="black")

        if annotations == True:
            ax[0].annotate(
                '{:.3f}'.format(fraction_np),
                (arg, fraction_np),
                size=12,
                ha='left' if kd_np >= kd_p else 'right',
                va='top' if kd_np >= kd_p else 'bottom',
                textcoords='offset points',
                xytext=(10 if kd_np >= kd_p else -10, 0),
                arrowprops={"arrowstyle": "- ", "shrinkA": 2,
                            "shrinkB": 2, "linewidth": 0.5},
            ).draggable()
            ax[0].annotate(
                '{:.3f}'.format(fraction_p),
                (arg, fraction_p),
                size=12,
                ha='right' if kd_np >= kd_p else 'left',
                va='bottom' if kd_np >= kd_p else 'top',
                textcoords='offset points',
                xytext=(-10 if kd_np >= kd_p else 10, 0),
                arrowprops={"arrowstyle": "- ", "shrinkA": 2,
                            "shrinkB": 2, "linewidth": 0.5},
            ).draggable()

    for i, arg in enumerate(args):
        ax_delta.vlines(arg, -1, 1, color='gray', linestyle=':')
        r = on_off_change(kd_np, kd_p, arg)
        fraction_diff = r['Fraction difference']
        fc = r['Fold-change']
        ax_delta.vlines(arg, 0, fraction_diff, color='gray', linestyle=':')
        globals()[f'a{n}'] = ax_delta.plot(
            [arg, arg],
            [0, fraction_diff],
            marker='o',
            ms=3.5,
            color="black",
            linestyle=''
        )
        n += 1
        ax_delta.plot(
            [arg, arg],
            [0, fraction_diff],
            color="black",
            linestyle='-'
        )
        globals()[f'a{n}'] = ax[1].plot(
            [arg], [fc], marker='o', ms=4, color='green')
        n += 1
        if annotations == True:
            if type == 'optimal':
                annotation = r'$\Delta_{max}$ = ' + \
                    '{:.3f}'.format(fraction_diff)
            else:
                annotation = f'$\Delta$ = '+'{:.3f}'.format(fraction_diff)
            ax[1].annotate(
                annotation,
                (arg, 10 ** (3*fraction_diff)),
                size=12,
                ha='center' if not type == 'optimal' else 'left',
                va='bottom' if kd_np >= kd_p else 'top',
                textcoords='offset points',
                xytext=(0, 10 if kd_np >= kd_p else -10),
                arrowprops={"arrowstyle": "- ", "shrinkA": 2,
                            "shrinkB": 2, "linewidth": 0.5},
                visible=False if type == 'baseline' else True,
            ).draggable()
            ax[1].annotate(
                'r = '+'{:.3f}'.format(fc),
                (arg, fc),
                size=12,
                xytext=(arg/10000, (0.05 * (i+1))
                        ) if kd_np >= kd_p else (arg/10000, (20 * (i+1))),
                arrowprops={"arrowstyle": "- ", "shrinkA": 2,
                            "shrinkB": 2, "linewidth": 0.5, "color": "green"}
            ).draggable()

        if type == 'optimal':

            r = kd_np/kd_p

            ax[1].annotate(
                r'r$_{max}$ = ' +
                '{:.3f}'.format(r) if annotations == True else None,
                (2e-14, r*1.4 if 500 > r >= 1 or r < 0.002 else r/2.5),
                size=12,
                ha='left'
            ).draggable()

    fig.tight_layout()
    if annotations == False:
        for i in range(1, n):
            mplcursors.cursor(globals()[f'a{i}'], multiple=True)
    plt.show()
