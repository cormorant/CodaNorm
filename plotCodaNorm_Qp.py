#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plot CodaNorm results for each Station.

Final version of resulting plot of CodaNorm method.
"""
import os
import sys
import numpy as np
from statsmodels.robust.robust_linear_model import RLM

import matplotlib.pyplot as plt
plt.style.use('classic')
#from qopen.util import gmean, gerr#gstat#robust_stat
from matplotlib import rcParams
FONTSIZE = 18
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = FONTSIZE


PATH = "values"

FREQS =  (1.5, 3., 6., 12., 24.)

LIMITS = (0.5, 1., 2., 4., 8.)


Vs = 3.51# km/s

MAX_DIST = 70


font = {
    'family': 'arial',
    'size': FONTSIZE,
}

bbox_dict = dict(
    facecolor='w',
    edgecolor="k",#'none',
    alpha=0.75
)

#DATE_E	TIME_E	LAT	LON	K	DIST	SNR	Ap	As	Ac	A100	A100RMS
DTYPE = [('DATE_E', "|U10"), ('TIME_E', "|U11"), 
    ('LAT', float), ('LON', float), ('K', float),
    ('DIST', float), ('SNR', float), 
    ('Ap', float), ('As', float), ('A100', float), ('A100RMS', float)]


def load_values(filename, station=None):
    """ load values from TXT-file with header:
    DATE_E TIME_E LAT LON K CHA ALPHA DIST REAL_DIST AZ SNR P_Z C_Z LN """
    
    data = np.loadtxt(filename, 
        delimiter="\t", skiprows=1, dtype=DTYPE)
    
    # hide EQs with dist > MAX_DIST
    ind = np.where( data['DIST'] <= MAX_DIST )
    data = data[ind]
    
    return data
    

def linear_fit(y, x, m=None, method='robust', Model=RLM):
    """Code from Qopen: Linear fit between x and y

    :param y,x: data
    :param m: fix slope at specific value
    :param method: one of ('least_squares', 'robust')
    :return: slope a and intercept b of y = ax + b
    """
    if m is None:
        X = np.empty((len(y), 2))
        X[:, 0] = x
        X[:, 1] = 1
        res = Model(y, X).fit()
        return res.params
    else:
        X = np.ones(len(y))
        res = Model(np.array(y) - m * np.array(x), X).fit()
        return m, res.params[0]


def calc_power_fit(freqs, values):
    """ calc Power fit of data
    But use log-log function, and linear fit instead """
    
    # take natural logarithm of both freqs and values
    freqs, values = np.log(freqs), np.log(values)
    
    # linear fit
    k, b = linear_fit(values, freqs)
    
    # k - frequency parameter `n`
    Q0 = np.exp(b)
    
    # output
    print('Q = %.2f * f ^ %.4f' % (Q0, k))
    
    # return Q0, n
    return Q0, k


def main(station, freq, ax):
    """ plot results from output files of CodaNorm_v2 """
    # NO! we have data for all channels in ONE file for every FREQuency
    FR = '%.2f' % freq
    
    # for evenry channel in N, E do
    # save Q values per channel
    Qs = []
    
    CHANNEL = "Z"
    path = os.path.join(PATH, station)#, CHANNEL)
    # find filename like: `TRGvalues_E_0.75Hz_SD10_T35.txt`
    
    files = [os.path.join(path, f) for f in os.listdir(path) if FR in f]
    # just files no dirs
    files = [f for f in files if os.path.isfile(f)]
    # must be 1 file for this channel and frequency
    if not len(files) == 1:
        print("Must be 1 file for FR %s" % freq)
        return 0, 0
    filename = files[0]
    
    # parse-read file
    try:
        data = load_values(filename, station=station)
    except IndexError as e:
        print(e)
        return
    
    #=== CALC values with Z(R) corr
    # make correction for Geom spreading...
    
    # distances array
    X = data['DIST']
    
    # geometrical spreading param for dist < 50 km is always == dist**-1
    Z_R = 1 / X
    # but for dist 50-70 may be 1/50 ???
    ind = np.where(X > 50)
    # make 1/50
    Z_R[ind] = 1/50
    
    #=======
    Y = np.log(data["Ap"] / ( data["A100RMS"] * Z_R ))
    #=======
    
    #===========================================================
    # ROBUST method Fit:
    a, b2 = linear_fit(Y, X)
    Y2 = a * X + b2
    
    # calc overall Q value
    # Q = Pi * f / (v * b)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Q = -np.pi * freq / (Vs * a)
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    Qs += [Q]
    
    # plot ROBUST results
    #===
    _label = "$Q_P = %.0f$" % Q
    
    Qline, = ax.plot(X, Y2, "-b", lw=1.5, 
        label=_label, zorder=222)
    #===
    
    # how much items
    _NUM = X.size
    print("CH = %s \t Freq = %.2f \t Q = %.0f \t N = %d" % (CHANNEL, freq, Q, _NUM))
    
    #===
    Qall, = ax.plot(X, Y, "ow", markersize=7,
        markeredgecolor="k", zorder=111) 
    
    ax.legend(loc='lower left', prop={"size":FONTSIZE})
    
    return Qs


if __name__ == '__main__':
    
    STATIONS = 'BGT FFNB KELR UZR'.split()
    
    #for STATION in STATIONS:
    STATION = "UZR"
    
    print("<<<<<")
    print(STATION)
    print(">>>>>")
    
    # save results here
    Qvalues = {}
    for _fr in FREQS: Qvalues[_fr] = []
    
    #=== Plot ln(As*R/Ac) vs. distance
    # prepare Figure
    fig, axes = plt.subplots(nrows=5, ncols=1, sharex=True, sharey=True, figsize=(9, 12), dpi=200)
    
    # Big title
    fig.suptitle(u'Station `%s`' % STATION, fontsize=FONTSIZE+2)
    
    (ax1, ax2, ax3, ax4, ax5) = axes
    
    #for ax in axes:
    ax3.set_ylabel(r"$ln(A_P \cdot Z(R)/A_C)$")
    
    ax5.set_xlabel(u"Distance, km")
    
    #=== do the job
    
    for FREQ, LIMIT, ax in zip(FREQS, LIMITS, axes):
        #print(FREQ, LIMIT)
        freq_min, freq_max = FREQ - LIMIT, FREQ + LIMIT
        ax.set_title(u"Frequency {0} Гц ({1}–{2} Гц)".format(("%.3f" % FREQ).rstrip('0').rstrip('.'), 
            ("%.3f" % freq_min).rstrip('0').rstrip('.'), ("%.3f" % freq_max).rstrip('0').rstrip('.')).replace('.0', ""))
        
        #===
        Qvalues2 = main(STATION, FREQ, ax)
        #===
        
        # save result Q for tis freq
        if Qvalues2 is not None:
            if len(Qvalues2)==1:
                Qvalues2 += [Qvalues2[0]]
            # add mean value
            Q1, Q2 = Qvalues2
            Qvalues[FREQ] += [(Q1+Q2)/2]
        
        # axis details
        #ax.set_xlim(20, 80)
        ax.set_xlim(0, 75)
        ax.set_ylim(0, 10)
    #===
    
    #plt.subplots_adjust(left=0.05, bottom=0.1, right=.98, top=0.9, wspace=0.05, hspace=0.15)
    plt.subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.97, 
        hspace=0.25, wspace=0.05)
    
    # save figure
    save_to_dir = os.path.join(PATH, STATION)
    if not os.path.exists(save_to_dir): os.makedirs(save_to_dir)
    outfilename = os.path.join(save_to_dir, "Qp_{}__70km_H.png".format(STATION))
    #plt.show()
    plt.savefig(outfilename)
    
    plt.close()
    
    #================================
    # Finally, make 1 resulting plot, for 20-40-60 overall Q value
    fig, ax = plt.subplots(figsize=(12, 9), dpi=200)

    # axis details
    ax.set_ylabel(r"$Q_P$")
    ax.set_xlabel(u"Частота, Гц")
    
    #=========================
    # do the job, collect and calc resulting Q value, for all SD (window length)
    
    Q_VALUES, Q_ERRORS_PLUS, Q_ERRORS_MINUS = [], [], []
    
    for FREQ, LIMIT in zip(FREQS, LIMITS):
        # values for this freq
        A = np.array(Qvalues[FREQ])
        # we need final 1 value with error (make errorbar!): yerr=ERROR
        QQ = int(np.round(A.mean()))
        
        Q_VALUES += [ QQ ]
        
        Q_ERRORS_PLUS += [A.max() - QQ]
        Q_ERRORS_MINUS += [np.abs(A.min() - QQ)]
        
        # maybe plot text values already?
        y_plus = A.max()-QQ
        y_plus2 = y_plus * 2.5
        y_plus2 = y_plus2+5 if FREQ==6. else y_plus2
        ax.text(FREQ, QQ+y_plus2, "$%d\pm%d$"%(QQ, y_plus), fontdict=font, bbox=bbox_dict, zorder=999)
    
    #=========================
    
    #=== plot Q values with error
    #First row contains the lower errors
    yerr = np.array([Q_ERRORS_MINUS, Q_ERRORS_PLUS]) # second row - upper
    
    # make power fit of Q values:
    Q0, nn = calc_power_fit(FREQS, Q_VALUES)
    # label of fit
    fit_label = '$Q = %.0f \cdot f ^ {%.3f}$' % (Q0, nn)
    
    ax.loglog(np.array(FREQS), Q_VALUES, "r:", lw=1.5, zorder=-1)
    
    # plot power-law fit
    _freqs = np.array(FREQS)
    Qfit = Q0 * np.power(_freqs, nn)
    ax.loglog(_freqs, Qfit, "-b", lw=1., label=fit_label)
    
    ax.errorbar(np.array(FREQS), Q_VALUES, 
        # Separate - and + values for each bar. 
        yerr=yerr, 
        fmt="wo", #marker="_", normally `o`
        markeredgecolor="k", ms=5,#markersize
        #uplims=True, lolims=True,
        ecolor='k', elinewidth=1, capsize=5,
    )
    
    ax.set_xlim(0.5, 50)
    ax.set_xticks(np.array(FREQS))
    ax.set_xticklabels([("%.2f"%_f).rstrip("0").rstrip(".") for _f in FREQS])
    
    ax.legend(loc='upper left', fancybox=True)#'lower right'
    
    plt.savefig('{0}/{1}/Q0_{1}_Qp_70km_H.png'.format(PATH, STATION))
    
    plt.close()

