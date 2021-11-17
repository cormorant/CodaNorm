#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division, print_function
#
#  CodaNorm.py
#  
#  Copyright 2014 petr <petr@SEISMOGRAMMA>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
"""
Program for calculating the RMS amplitudes of bulk waves and coda waves. 
Input data format  - MSEED.

UPDATE 2020-06-25:
    - added command line arguments, in order to override "coda.conf" params

TODO: use json for config, also for results saving
"""
APP_NAME = "CodaNorm"
__version__="0.3"
COMPANY_NAME = 'GIN SB RAS'

import os
import sys
import time
import datetime
import math

#import obspy
from obspy import UTCDateTime, Trace, Stream, read as read_mseed
from obspy.geodetics.base import gps2dist_azimuth

import numpy as np
import scipy

# CodaNormlib code
from CodaNormLib import (STATIONS, linear_fit_qopen, read_settings, 
    RMS, calc_signal_noise_ratio, get_waveforms, get_value_for_window)

#import json
# command line arguments may override default settings from `coda.conf` file
import argparse

# grapihcs
import matplotlib.pyplot as plt
plt.style.use('classic')
from matplotlib import rcParams
FONTSIZE = 16
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = FONTSIZE-2


def module_path():
    ''' add link to where this coda is stolen '''
    if hasattr(sys, "frozen"):
        return os.path.dirname(sys.executable)
    return os.path.dirname(__file__)

# get current dir, may vary if run from EXE file
CurrDir = module_path()

#sys.stdout = open(os.path.join(CurrDir, "log", "CodaNorm.log"), "a")
#sys.stderr = open(os.path.join(CurrDir, "log", "CodaNorm.err"), "a")

CHANNELS = "N E Z".split()

CONFIG_FILENAME = os.path.join(CurrDir, "coda.conf")

#===============================================================================


def calculate(Freq, f1, f2, stream, 
    seconds_E, seconds_P, seconds_S, 
    dt, dt_P, dt_S, Settings, dist, Acoda=None):#azBA=0):
    """ 
    Freq, f1, f2
    stream - a copy of 'slices' original STREAM object that saves traces
    
    seconds_E - seconds relative to the beginning of the file - time at source
    (should always == 0 if REQUEST_WINDOW = (0, ...))
    
    seconds_S - seconds S-wave from the T0
    
    DT, DT_P, DT_S - UTCDateTime objects
    Settings 
    dist - epi distance, km
    Acoda - A100 or A40 etc
    """
    
    # save original trace
    original_trace = stream[0].copy()
    T0 = original_trace.stats.starttime
    
    sr = float("%.0f" % original_trace.stats.sampling_rate) # 100 or 20 Hz normally
    
    if sr == 50.:
        if Settings["verbose"]:
            print("\nLets not use Irkut-24 digitizer data")
            return
    
    #===
    # according to [Aki 1980] [Pavlenko 2008] [etc.]
    # нужно брать коду за фиксированныое время - например 40 с после Т0 (время в очаге)
    # однако для расстояний > 70 км, 40 сек будет раньше времени 2*Ts !
    
    #===
    # time of Coda start
    dt_Coda = dt + Acoda
    
    # how much is coda start and end in seconds
    # cut Coda from Coda1 to Coda2 - -10s back -- +10s forward (if SD == 20)
    Coda1, Coda2 = Acoda - SD/2, Acoda + SD/2
    
    # time of coda start 
    dt_Coda1 = dt + Coda1
    # and End
    dt_Coda2 = dt + Coda2
    
    # check, maybe A100 < 2*Ts, then skip it
    if Acoda < (seconds_S * 2):
        if Settings["verbose"]:
            print("\nA`%d` amplitude is too earlier then `%.0f`!" % (Acoda, Coda1))
        # exit only if diff > 10-15 s, cause we can use anyway 1.5*Ts
        if Acoda < (seconds_S * 2 + 10):
            return
    
    # calc SD1 manually (if not set)
    # or use fixed length for P- and S-waves
    SD1 = Settings['sd1']
    if SD1 == 0:
        # window for Direct body-waves = 0.25 * Ts
        SD1 = int( math.ceil(0.25 * seconds_S) )

    # P-wave window is half of S-window
    SD1P = int(math.ceil(SD1 / 2))
    
    # check: S-window and Coda-window may intersects!
    if (Coda1 - seconds_S) < SD1:
        if Settings["verbose"]:
            print("\nWindow S and Coda intersects for event `%s`!" % dt)
        return
    
    #===
    # save resulting values here
    result = []
    
    # filtering
    filter_params = {'freqmin': f1, 'freqmax': f2, 'type': 'bandpass', 
        'corners': Settings["corners"], 'zerophase': True}
    # Butterworth filter in place
    stream.filter(**filter_params)
    
    if Settings["component"] == "H":
        # get ?N and ?E channels
        __N = stream.select(channel='?HN')[0].data
        __E = stream.select(channel='?HE')[0].data
        # its Hilbert transform
        hilbN = scipy.fftpack.hilbert(__N)
        hilbE = scipy.fftpack.hilbert(__E)
        # calc resulting components
        __N = np.power(__N, 2) + np.power(hilbN, 2)
        __N = np.power(__N, 0.5)
        
        __E = np.power(__E, 2) + np.power(hilbE, 2)
        __E = np.power(__E, 0.5)
        
        # and it's Average over two channels (result of calculations)
        RMS_component = np.mean( np.array([ __N, __E ]), axis=0 )
        # remove 2 components
        for _ in range(2): stream.pop(0)
        stream[0].data = RMS_component
    elif Settings["component"] in ("N", "E", "Z"):# use just 1 channel
        # channel data
        __Z = stream.select(channel='?HZ')[0].data
        # its Hilbert transform
        hilbZ = scipy.fftpack.hilbert(__Z)
        # resulting envelope
        __Z = np.power(__Z, 2) + np.power(hilbZ, 2)
        # square root from it
        stream[0].data = np.power(__Z, 0.5)
    stream[0].stats.channel = "RMS"
    
    #=== start calculating
    # for Qp we use only Z-channel!!! for Qs - selected...
    assert stream.count() == 1, "\nMust be 1 trace! we have: %s" % stream
    
    trace = stream[0]
    
    # check SNR on filtered signal!
    SNR, RMS_NOISE = calc_signal_noise_ratio(trace, sr, dt, dt_Coda2, 3)# 3 sec
    
    # save SNR
    result += [SNR]
    
    if SNR < Settings['minsnr']:
        if Settings["verbose"]:
            print("\nSNR for freq %s is too SMALL (%.1f)!..." % (Freq, SNR))
        return result
    
    #===================================================================
    #=== Calculations
    
    # get envelope around 40 (or 100) second, aorund Acoda
    envelope_part = trace.slice(starttime=dt_Coda1, endtime=dt_Coda2)
    envelope_times = envelope_part.times() + Coda1
    
    # calc A100 value
    # just mean value like Eulenfeld does: tr.data = tr.data / np.mean(sl.data)
    y100 = np.mean(envelope_part.data) # mean value
    y100_RMS = RMS(envelope_part.data) # RMS
    
    #===
    # get filtered parts of data
    
    # P-window data
    # check that P-window and S-window do not intersects!
    if (seconds_S - seconds_P) <= SD1P:
        if Settings["verbose"]:
            print("\nWindow S and P intersects for %s!" % dt)
        Pwindow = None
    else:
        Pwindow = stream.slice(dt_P, dt_P+SD1P)[0].data
    
    # S-window data
    Swindow = stream.slice(dt_S, dt_S+SD1)[0].data
    
    #===========================================================================
    # *** Crucial place for calculations ***
    
    # result for P-bulk-window
    P_window_result = get_value_for_window(Pwindow, maximum=Settings["max"])
    # S-window
    S_window_result = get_value_for_window(Swindow, maximum=Settings["max"])
    
    #TODO: remove coda_value, we calc RMS(Acoda), A100, A40 etc
    coda_value = 0#get_value_for_window(coda_window, maximum=Settings["max"])
    
    # save results +A100 +slope koeff
    result += [P_window_result, S_window_result, coda_value, y100, y100_RMS]
    
    #===========================================================================
    #=== calculations Ended
    
    
    #===========================================================================
    # Move all plotting issues here, move out of calculating code
    #===========================================================================
    
    if PLOT:
        # make Figure
        fig, (ax1, ax2) = plt.subplots(figsize=(12, 7), nrows=2, ncols=1, sharex=True)#, dpi=200)
        
        EVENT_NAME = "Dist = {:.0f} km | Stream {}--{} | Freq {} Hz ({}-{})".format(dist, 
            stream[0].stats.starttime, stream[0].stats.endtime, Freq, f1, f2)
        fig.suptitle(EVENT_NAME, fontsize=FONTSIZE)
        
        # plot original signal
        ax1.plot(original_trace.times(), np.abs(original_trace.data), 
            color="grey", lw=.5, zorder=111, label="NT.%s.00.BHN+BHE" % STATION)
        # plot filtered signal
        ax1.plot(trace.times(), trace.data, "b", lw=1., zorder=222, label="RMS-envelope")
        
        # plot filtered
        #labelFiltered = 'Фильтр {}-{} Гц'.format(f1, f2)
        ax2.semilogy(trace.times(), trace.data, "k", lw=0.5, zorder=222, alpha=0.75)#label=labelFiltered, 
        
        # mark time of Event, P and S time by vertical lines, start -- end
        for ax in (ax1, ax2):
            ax.axvline(x=seconds_E, linestyle="--", color="y", lw=2.) # Event
            ax.axvline(x=seconds_P, linestyle="--", color="b", lw=0.5) # P
            ax.axvline(x=seconds_P+SD1P, linestyle="--", color="b", lw=0.5) # P+SD1p
            # S-window
            ax.axvline(x=seconds_S, linestyle="--", color="k", lw=1) # S
            ax.axvline(x=seconds_S+SD1, linestyle="--", color="k", lw=1) # S+SD1
            # mark coda
            ax.axvline(x=Coda1, linestyle="--", color="r", lw=2.) # coda
            ax.axvline(x=Coda2, linestyle="--", color="r", lw=2) # coda
            # mark noise RMS_NOISE
            ax.axhline(y=RMS_NOISE, linestyle="--", color="c", lw=.5) # noiz
        
        ax1.plot(Acoda, y100, "om", markersize=7, markeredgecolor="k", zorder=9999)
        ax2.plot(Acoda, y100, "om", markersize=7, markeredgecolor="k", zorder=9999)
        
        # mark result for S-window
        ax2.plot(seconds_S, S_window_result, "om", markersize=7, 
            markeredgecolor="k", zorder=999)
        
        # plot ec. details
        y_pos = trace.data.max() / 2
        
        #=== selected windows
        # P
        P_time = trace.times()[:Pwindow.size] + seconds_P
        ax2.plot(P_time, Pwindow, "-r", lw=1.)
        
        # plot S-window (Swindow)
        S_time = trace.times()[:Swindow.size] + seconds_S
        ax2.plot(S_time, Swindow, "-r", lw=1.)
        
        # Coda (envelope) on second axis
        ax2.plot(envelope_times, envelope_part.data, "-r", lw=1.)
        
        # plot SNR value text
        ax1.text(x=seconds_E+0.5, y=y_pos, s="%.1f" % SNR) # SNR
        
        # maximum Coda amplitude
        ax2.text(x=Coda1, y=y100_RMS, s="%.1f" % y100_RMS) #coda MAX
        
        # axis limits
        ax2.set_xlim(seconds_E - 5, Coda2 + 5)
        # y
        ax1.set_ylim(-5, original_trace.data.max())
        ax2.set_ylim(None, trace.data.max()*1.5)
        
        # legend
        ax1.legend(loc="best")
        #ax2.legend(loc="best")#upper right
        
        plt.tight_layout()
    
        if Settings["show"]:
            plt.show()
        else:
            # just save Figure
            save_to_dir = os.path.join("img", STATION, str(stream[0].stats.starttime)[:19].replace(":", "-"))
            if not os.path.exists(save_to_dir): os.makedirs(save_to_dir)
            outfilename = os.path.join(save_to_dir, "{}_{}_{}__{}s.png".format(STATION, Settings["component"], Freq, SD))
            plt.savefig(outfilename)

    # the end
    return result


def parse_command_line_args(args=None):
    """ parse command line arguments """
    description = 'CodaNorm arguments that override default settings'
    parser = argparse.ArgumentParser(description=description)
    
    # argument on which channel to calc
    parser.add_argument("-c", "--channel", help="Channel to use (N, E, Z).", 
        choices=CHANNELS, nargs='+')
    
    # add option to calc with another coda window length
    parser.add_argument("--sd", type=int, default=0,
        help="Specify length of Coda-window (in seconds)")
    
    # maybe we want to plot anyway
    parser.add_argument('--plot', help='Plot results', action='store_true')
    
    # change station
    parser.add_argument("--station", default=None,
        help="Specify station to use")
    
    # verbose mode
    parser.add_argument("-v", "--verbose", help="output verbosity", 
        action="store_true")
    
    args = parser.parse_args(args)
    
    return args


def get_geometrical_spreading_factor(dist, R1=60, R2=100, power2=0.5):
    """ 3-segment function of Geometrical spreading factor:
    0 - R1  -- 1 / R
    R1 - R2 -- 1 / R1
    > R2    -- R ** -0.5 """
    
    # for distances <= 1.5 * Moho, direct waves adn geom spreading factor == -1
    if dist <= R1:
        Z = 1 / dist
    elif R1 < dist <= R2:
        Z = 1 / R1
    else:
        # > R2 (100 km)
        Z = 1 / R1 * (dist / R2)**-power2
    
    return Z


if __name__ == '__main__':
    #===
    # Для каждой записи получить:
    # времена события (задается во входном каталоге)
    # - время вступления волны P для данной станции
    # - время вступления волны S для данной станции
    # - имя файла для обработки, откуда считывать данные
    #=== считывание файла настроек
    try:
        Settings = read_settings(CONFIG_FILENAME)
    except BaseException as e:
        print(e)
        sys.exit(1)
    
    #=== Parse settings/params
    
    # === may override Settings by command line arguments
    # parse command line arguments
    args = parse_command_line_args()
    
    # change Plot or not
    if args.plot: Settings["plot"] = args.plot
    # change channel if needed
    if args.channel is not None:
        Settings['component'] = args.channel[0]
    
    # may want another SD:
    if (args.sd is not None) and (args.sd != 0):
        Settings["sd"] = args.sd
        SD = args.sd
    
    # verbose or not
    Settings["verbose"] = True if args.verbose else False
    
    # ======
    
    # get name if input file
    STATION = Settings["station"] if args.station is None else args.station
    if args.station is not None: Settings["station"] = args.station

    # get stations coords, Input filename and Data MSEED volume
    Settings['station_lat'], Settings['station_lon'], InputFile, DATAFILE = STATIONS[STATION]

    # global variables here
    Freqs, Limits = Settings["freqs"][::2], Settings["freqs"][1::2]

    # etc Settings
    PLOT = Settings["plot"]

    # total length of Coda-window (seconds)
    SD = Settings["sd"] # LENGTH OF CODA WINDOW

    # total length of Body-wave bulk window (seconds)
    SD1 = Settings["sd1"]

    # waves velocities and it's relation
    Vp = Settings["vp"]
    Vs = Settings["vs"]
    VpVs = Vp / Vs

    # what data to fetch (we will use only 1 channel)
    if Settings["component"] == "H":
        which_channel = "?H?"
    elif Settings["component"] in CHANNELS:
        which_channel = "?H%s" % Settings["component"]# if not Settings["rms"] else "?")

    # NETWORK station etc
    net, sta, loc, cha = "??", STATION, "00", which_channel

    # read data Stream
    if Settings["verbose"]:
        print("Read MSEED data from `%s`..." % DATAFILE)
    DATAFILE = os.path.join("data", STATION, DATAFILE)
    # will cache overall data in memory? read overall data onto memory
    STREAM = read_mseed(DATAFILE, format="MSEED")
    
    # input file in dir "catalog/[STA]/"
    InputFile = os.path.join("catalogs", STATION, InputFile)
    if not os.path.exists(InputFile):
        print("Input file '%s' not found!" % InputFile)
    
    if Settings["verbose"]:
        print("Read events from `%s`..." % InputFile)
    # load input file
    Catalog = np.genfromtxt(InputFile, dtype=None, autostrip=True, names=True,
        loose=True, invalid_raise=False, encoding="utf-8")
    # имена столбцов
    names = Catalog.dtype.names
    # проверить что все нужные столбцы в файле есть
    for col in ("DATE_E", "TIME_E", "TIME_P", "TIME_S"):
        if not col in names:
            print("Coulmn '{}' not found in input file. Fix it!".format(col))
            sys.exit(1)
    # columns: DIST_KM AZIMUTH
    names = list(names) + ["DIST_KM", "SNR", "AZIMUTH_AB", "AZIMUTH_BA"]
    
    # ready to start, show settings
    
    if Settings["verbose"]:
        print("\tSettings:\n")
        for k, v in Settings.items(): print(k, v)
        
        print("\tCommand line arguments:")
        print(args)
        
        # nice wait
        _N = 3
        while True:
            if _N <= 0: break
            # nice output
            s = "Wait %d seconds..." % _N
            sys.stdout.write("\r" + s)
            sys.stdout.flush()
            #
            time.sleep(1)
            _N -= 1
        print()
        print("\t...started")
    
    #===========================================================================
    
    print("*"*33)
    print("\t Starting program `{}` v.{}. Station data `{}`".format(APP_NAME, __version__, STATION))
    
    # start iterating over Catalog
    for nnn, line in enumerate(Catalog):
        # nice output
        s = "%s of %s\t%s" % (nnn+1, Catalog.size, " ".join(map(str, line))[:37])
        sys.stdout.write("\r" + s)
        sys.stdout.flush()
        
        # обрабатывать данные по событию
        parts = dict(zip(names, line))
        
        # epicenter coords
        lat, lon = parts["LAT"], parts["LON"]
        
        # calc distance and azimuthes from epicenter
        try:
            dist, azAB, azBA = gps2dist_azimuth(Settings["station_lat"],
                Settings["station_lon"], lat, lon)
        except ValueError as e:
            print(e)
            dist = 9999
        # Distance must be in kilometers (not in m)
        dist /= 1000.
        
        # hypocenter correction, use 15km default
        dist = (dist**2 + 15**2)**0.5
        
        # Distance checking...
        if not (Settings["mindist"] <= dist <= Settings["maxdist"]):
            if Settings["verbose"]:
                print("\n\tSkip event at epi distance %.0f km" % dist)
            continue
        
        # will recalc Coda amplitude to this abolute time (for example A100)
        Acoda = Settings["Acoda"]
        
        #== make UTC datetime(s)
        # time in source
        DT = UTCDateTime("%sT%s" % (parts["DATE_E"], parts["TIME_E"]))
        
        # time-P maybe `None` somewhere
        if parts["TIME_P"].upper() == "NONE":
            if Settings["verbose"]: print("\nERROR: Must have P-onset")
            #continue
        else:
            DT_P = UTCDateTime("%sT%s" % (parts["DATE_E"], parts["TIME_P"]))
            # P-wave onset in seconds
            seconds_P = DT_P - DT
        
        # may be no value for S-onset
        if parts["TIME_S"] == "None":# and TIME_P is not NONE
            # calc from time-P
            DT_S = DT + (DT_P - DT) * VpVs
        else:
            DT_S = UTCDateTime("%sT%s" % (parts["DATE_E"], parts["TIME_S"]))
            # S-wave onset in seconds
            seconds_S = DT_S - DT
            # also may have TIME_P == NONE
            if parts['TIME_P'].upper() == "NONE":
                # calc Tp from Ts
                seconds_P = seconds_S - dist / 8 # from formula (Ts-Tp)*8=dist
                DT_P = DT + seconds_P
        
        # check time valid
        if not (DT <= DT_P <= DT_S):
            print("\n Error (%s <= %s <= %s) 'DT <= DT_P <= DT_S' for %s" % (DT, DT_P, DT_S, line))
            continue
        
        #=== P-, S-onset from T0 in seconds?
        # time in source - 0
        seconds_E = 0
        
        # check Vp/Vs
        # from given time, calc velocity again
        calculatedVpVs = (dist/seconds_P) / (dist/seconds_S)
        
        if not (1.5 <= calculatedVpVs <= 2.):
            if Settings["verbose"]:
                msg = "Vp/Vs must be in range 1.5--2. (Vp/Vs=`%.2f`)" % calculatedVpVs
                print(msg)
                continue
        
        #=== Fetch data: request this waveforms from MSEED volume (in seconds)
        
        # calc seismogram length needed
        SEISMOGRAMM_MAX_LEN = Acoda + SD + 10
        
        # we will fetch data for this times
        REQUEST_WINDOW = (0, SEISMOGRAMM_MAX_LEN)
        
        # Retrieve data
        t1, t2 = DT + REQUEST_WINDOW[0], DT + REQUEST_WINDOW[1]
        
        kwargs2 = {"stream": STREAM, 'network': net, 'station': sta, 'location': loc,
            'channel': cha, 'starttime': t1, 'endtime': t2, 'event': 'event'}
        
        # Fetch data
        stream = get_waveforms(**kwargs2)
        
        if stream.count() == 0:
            if Settings["verbose"]: print("\nNo data for event `%s`" % DT)
            continue
        
        # maybe duplicate files?
        if stream.count() == 6: stream.merge(method=1)
        
        # must be all 3 components
        if Settings['component'] == "H":
            if not stream.count() == 3:
                print("\nERROR: must be 3 channel to make RMS compomnent!")
                continue
        else:
            err_msg = "\nmust be 1 ch.: %s (%s)" % (stream, which_channel)
            if stream.count() != 1:
                print(err_msg)
                continue
        
        # detrend (this includes removing mean value)
        stream.detrend("linear")
        
        #=============
        #== считать окна, для времен вступлений, 1 раз - для всех частот
        # начало и конец файла
        T0, Tend = stream[0].stats.starttime, stream[0].stats.endtime
        
        # must check that:
        # 1) P- and S-window is here: Time of P- and S- start and end we have
        # 2) Coda window is all in this data: start and end
        
        # DT DT_P DT_S - UTCatetime object that saves start of this waves (DT id T0)
        
        # check that P-wave window we have here...
        if (DT_P < T0) or (DT_P > Tend):
            # if one of this - not enough data!!!
            err_msg = "\nP-window must start after T0 and end before Tend"
            if Settings["verbose"]: print(err_msg)
            continue
        
        # check that S-window we have overall
        if (DT_S < T0) or (DT_S > Tend):
            # if one of this - not enough data!!!
            err_msg = "\nS-window must starts after `T0` and ends before `Tend`"
            if Settings["verbose"]: print(err_msg)
            continue
        
        # Main check: Coda-window 
        # start of Coda-window must be also after T0 and before Tend
        start_of_Coda = DT + Acoda
        
        if (start_of_Coda < T0) or (start_of_Coda > Tend):
            err_mgs = "\nCoda-window start must be between `T0` and `Tend`"
            if Settings["verbose"]: print()
            continue
        # check end of Coda
        end_of_Coda_window = start_of_Coda + SD
        if (end_of_Coda_window > Tend):
            if Settings["verbose"]: print("\nCoda-window start must endbefore `Tend`!")
            continue
        
        # check we will make: T0 must be EXACTLY what we fetched...
        if not T0 == t1:
            # first check the difference...
            diff = abs(T0 - t1)
            # if difference < 0.1 second, ok
            if diff < 0.1: pass
            else:
                if Settings["verbose"]: 
                    print("\nStream must start at the same time as T1 (%s != %s)! \nDiff = %s" % (T0, t1, diff))
                continue
        
        #===
        # WARNING! recalc distance according to formula 8 * (Ts - Tp + 0.2)
        #real_dist = 8 * (seconds_S - seconds_P + 0.2)
        #if Settings["verbose"]:
        #    print("\nReal_dist = %.2f, dist was %.2f" % (real_dist, dist))
        #dist = real_dist
        
        #===============================================================
        # calculation for all frequencies
        
        calc_freq_values = []
        
        for NumFreq, (Freq, plus_min) in enumerate(zip( Freqs, Limits )):
            # подсчитать границы от центральной частоты, low and high freq corners
            f1, f2 = Freq - plus_min, Freq + plus_min

            #===========================================================
            try:
                # функция "calculate" возвращает 3 значения:
                # 0 - SNR
                # 1 - value (RMS or maxima) in body-wave window
                # 2 - value (RMS or maxima) in Coda-wave window
                result = calculate(Freq, f1, f2, stream.copy(),
                    seconds_E, seconds_P, seconds_S, 
                    DT, DT_P, DT_S,
                    Settings, dist,#if real_dist == 0 will use dist
                    Acoda
                )
                plt.close()
            except KeyboardInterrupt:
                print("\nInterrupted by user...")
                PLOT = False
                #break
            
            #===========================================================
            # if no values return, skip this freq...
            if (result is None) or ( len(result) == 1 ):
                plt.close()
                # надо считать по другим частотам, если по 1-й неудачно? Yeas!
                calc_freq_values.append( 0 )
                continue
            
            # parse results
            # SNR, P-value, S-value, Coda-value, A100, slope
            SNR, P_value, S_value, C_value, Y100, Y100RMS = result
            
            # save distance-Amplitude - Ap/As/Ac/Y100 values
            datapath = os.path.join("values", STATION)
            if not os.path.exists(datapath): os.makedirs(datapath)
            # filename: STA_FREQ_SD_vals.txt
            datafilename = "{0}values_{1:0.2f}_SD{2}_T{3}.txt".format(
                STATION, Freq, SD, 0 if Acoda is None else Acoda)
            datafilename = os.path.join(datapath, datafilename)
            # Do not write all columns, just DIST and Amplitudes are necessary
            if not os.path.exists(datafilename):
                # create datafile and add Header(s)
                _dataf = open(datafilename, "w")
                _headers = "\t".join("DATE_E TIME_E LAT LON K DIST SNR Ap As Ac A100 A100RMS".split())
                _dataf.write(_headers + "\n")
                _dataf.close()
            
            with open(datafilename, "a") as _dataf:
                # write idE, idDir, DATE_E, TIME_E, LAT, LON, K, CHA, SNR, P_value, S_value_ C_value, C_Z, LN
                s = "{DATE_E}\t{TIME_E}\t{LAT:.3f}\t{LON:.3f}\t{K:.1f}\t".format(**parts)
                
                # distance: real and for calc. + Azimuth
                s += "%.2f\t" % dist
                
                # append results (SNR, P_value, S_value, C_value, A100) _k
                s += "\t".join(["%.5f" % _r for _r in result])
                
                _dataf.write(s + "\n")
                _dataf.close()

            #=== calculations ended
            #===========================================================
    
    # finished for all events from Catalog
    #===========================================================================
    
    print("\nDone!")
    print(":"*33)
