#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
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

APP_NAME = "CodaNorm"
__version__="0.2.0.2"
COMPANY_NAME = 'GIN SB RAS'

import os
import sys
import time
import datetime
import math
import warnings
warnings.filterwarnings("ignore", message='Comparing UTCDateTime objects of different precision is not defined will raise an Exception in a future version of obspy')

import obspy
from obspy import UTCDateTime, Trace, Stream, read as read_mseed
from obspy.geodetics.base import gps2dist_azimuth

import numpy as np

from scipy.signal import find_peaks
from scipy.stats.stats import pearsonr

# CodaNormlib code
from CodaNormLib import (STATIONS, linear_fit, read_settings, RMS, 
    calc_signal_noise_ratio, get_waveforms, remove_response)

# command line arguments may override default settings from `coda.conf` file
import argparse

# grapihcs
import matplotlib.pyplot as plt


CHANNELS = ("N", "E", "Z")

CONST_Tp_Ts = 1.76 # Vp/Vs constant for region, about 1.7-1.8


def module_path():
    ''' add link to where this coda is stolen '''
    if hasattr(sys, "frozen"):
        return os.path.dirname(sys.executable)
    return os.path.dirname(__file__)


# get current dir, may vary if run from EXE file
CurrDir = module_path()

CONFIG_FILENAME = os.path.join(CurrDir, "coda.conf")

#===============================================================================
def make_approximation(freq, tPeak, Peaks, T1, T2, t100, correct_to_Max=False):
    _y = np.log(Peaks)
    try:
        p = np.polyfit(tPeak, _y, 1)# we used natural logarithm
    except TypeError:
        return
    
    # what we did: p[0] - наклон, p[1] - вес
    k, b = p
    
    # y ≈ exp(-b) * exp(k * x)
    new_x = np.arange(T1, T2+1)
    new_y = np.exp(b) * np.exp(k * new_x)
    
    # also approximate from 0 to new_x[0]
    new_x0 = np.arange(0, new_x[0])
    new_y0 = np.exp(b) * np.exp(k * new_x0)
    
    # also calc A100
    y100 = np.exp(b) * np.exp(k * t100)
    
    #=== fix Y-values, to make max(Peaks) == max(Y-value)
    if correct_to_Max:
        # find mean of first 5-maximums, not to use maybe 1 outliner
        #max_Peak = Peaks.max()
        ind = np.argsort(Peaks)
        max_Peak = np.median(Peaks[ind][-5:])
        
        max_Y = new_y.max()
        ###if max_Peak > max_Y:
        diff = max_Peak - max_Y
        new_y += diff
        new_y0+= diff
        # also += Y100 value
        y100 += diff
        
    return new_x, new_y, y100, k, b, new_x0, new_y0


def calculate(Freq, f1, f2, stream, 
    seconds_E, seconds_P, seconds_S, 
    dt, dt_P, dt_S, Settings, dist, Tcoda=None):#azBA=0):
    """ 
    Freq, f1, f2
    stream - a copy of 'slices' original STREAM object that saves all data
    seconds_E - сколько это в секундах относительно начала файла - время в очаге
    (всегда 0 если у нас REQUEST_WINDOW = (0, ...))
    
    seconds_S - сколько секунд до S-волны от начала файла (или T0) """
    
    # save original trace
    original_trace = stream[0].copy()
    T0 = original_trace.stats.starttime
    
    sr = float("%.1f" % original_trace.stats.sampling_rate) # 100 Hz normally
    
    if sr == 50.:
        if Settings["verbose"]:
            print("Lets not use Irkut-24 digitizer data!!")
            return
    
    # stream must start at T0
    if abs(seconds_E - 0.) > 1:
        assert seconds_E == 0. , "\nseconds_E '%s' !=5" % (seconds_E)
    
    #===
    # according to [Aki 1980] [Pavlenko 2008] [Indus etc.]
    # нужно брать коду за фиксированныое время - например 40 с после Т0 (время в очаге)
    # однако для расстояний > 70 км, 40 сек будет раньше времени 2*Ts !!!
    
    # time of Coda
    dt_Coda = dt + Settings["koef"] * (dt_S - dt)
    
    # seconds (начало и конец коды)
    seconds_Coda1 = dt_Coda - T0
    seconds_Coda2 = seconds_Coda1 + SD
    
    # calc SD1 manually (if not set)
    SD1 = Settings['sd1']
    if SD1 == 0:
        # window for Direct body-waves = 0.25 * Ts
        SD1 = int( round(0.25*seconds_S) )
        
        SD1P = int(round(SD1 / 2))
    
    # check: S-window and Coda-window may intersects!
    if (seconds_Coda1 - seconds_S) < SD1:
        if Settings["verbose"]:
            print("\nWindow S and Coda intersects for event `%s`!!!" % dt)
        return
    
    #============
    # make Figure
    fig, (ax1, ax) = plt.subplots(figsize=(21, 9), nrows=2, ncols=1, sharex=True)
    
    EVENT_NAME = "Dist = {:.0f} | Stream {}--{}. Freq {} Hz ({}-{})".format(dist, 
        stream[0].stats.starttime, stream[0].stats.endtime, Freq, f1, f2)
    fig.suptitle(EVENT_NAME, fontsize=14)

    # будем считать
    result = []
    
    # filtering
    filter_params = {'freqmin': f1, 'freqmax': f2, 'type': 'bandpass', 
        'corners': Settings["corners"], 'zerophase': True}
    # Butterworth filter in place
    stream.filter(**filter_params)
    
    #=== start calculating
    # for Qp we use only Z-channel!!! for Qs - selected...
    assert stream.count() == 1
    
    trace = stream[0]
    
    #===
    # plot original signal
    ax1.plot(original_trace.times(), original_trace.data, "g", lw=1., zorder=111,
        label=trace.stats.channel)
    # plot signal before filtering
    ax1.plot(trace.times(), trace.data, "k", lw=0.5, zorder=222)#, alpha=0.75)
    
    # will use absolute values
    trace.data = np.abs(trace.data)

    # useless check
    if trace.data.max() == trace.data.min():
        print("\nWarning: may we have here solid lines!!!!!!!!")
        return
    
    # check SNR on filtered signal!
    # noise windows will always be = 3 second
    # calc at the end of the Coda-window!
    #TODO: for small Epi-distances (less 20-30 km), must use smaller windows!
    #TODO: must use several windows and see Noise at all this 
    SNR = calc_signal_noise_ratio(trace, sr, dt, dt_Coda+SD, 3)
    
    # save SNR
    result += [SNR]
    
    if SNR < Settings['minsnr']:
        if Settings["verbose"]:
            print("\nSNR for freq %s is too SMALL (%.1f)!..." % (Freq, SNR))
        return result
    
    # plot filtered
    ax.plot(trace.times(), trace.data, "k", lw=1, zorder=222)
    
    # mark time of Event, P and S time by vertical lines, start -- end
    for _ax in [ax1, ax]:
        _ax.axvline(x=seconds_E, linestyle="--", color="y", lw=2.) # Event
        _ax.axvline(x=seconds_P, linestyle="--", color="b", lw=0.5) # P
        _ax.axvline(x=seconds_P+SD1P, linestyle="--", color="b", lw=0.5) # P+SD1p
        # S-window
        _ax.axvline(x=seconds_S, linestyle="--", color="k", lw=1) # S
        _ax.axvline(x=seconds_S+SD1, linestyle="--", color="k", lw=1) # S+SD1
        # mark coda
        _ax.axvline(x=seconds_Coda1, linestyle="-", color="r", lw=2.) # coda
        _ax.axvline(x=seconds_Coda2, linestyle="-", color="r", lw=2) # coda
    
    
    #===================================================================
    #=== все вычисления по данной компоненте
    
    #=== Calculations here
    
    # get coda window to make envelope
    envelope_part = trace.slice(starttime=dt_Coda, endtime=dt_Coda+SD)
    envelope_times = envelope_part.times() + seconds_Coda1
    
    # get Peaks
    window_len = int( math.ceil((2. / Freq * sr)/10) * 10 ) + 1# odd number
    ind_peaks, _ = find_peaks(envelope_part, distance=window_len)
    
    tPeak, Peaks = envelope_times[ind_peaks], envelope_part[ind_peaks]
    # fin MAXIMUM from peaks
    maxPeak = Peaks.max()
    
    #===
    # make linear approximation of logarithm values
    result_approx = make_approximation(Freq, tPeak, Peaks, seconds_Coda1, seconds_Coda2, Tcoda)
    if result_approx is None:
        print("Error approximation for Freq `%.3f` and Event `%s`" % (Freq, EVENT_NAME))
    else:
        new_x, new_y, y100, _k, _b, new_x0, new_y0 = result_approx
        
    # plot approximation
    _label = r'$y=%.5f * exp(%.5f * x)$' % (np.exp(_b), _k)
    
    ax.plot(tPeak, Peaks, "*r", markersize=5, alpha=0.75, zorder=99999)
    ax.plot(new_x, new_y, "-r", lw=2., zorder=7777, label=_label)
    # result of approximation from 0 to new_x[0]
    ax.plot(new_x0, new_y0, ":r", lw=1.2, zorder=7777)
    
    
    ax.plot(Tcoda, y100, "or", markersize=12, markeredgecolor="w", alpha=0.5, zorder=9999)
    ax.legend(loc="upper right")
    ax1.legend(loc="upper right")
    
    #===
    # получить сами массивы данных (отфильтрованных)
    
    # P-window data
    # check that P-window and S-window do not intersects!
    if (seconds_S - seconds_P) <= SD1P:
        #if Settings["verbose"]:
        print("\nWindow S and P intersects for %s!!!" % dt)
        #return
        Pwindow = None
    else:
        Pwindow = stream.slice(dt_P, dt_P+SD1P)[0].data
    
    # S-window data
    Swindow = stream.slice(dt_S, dt_S+SD1)[0].data
    
    # Coda-window data
    try:
        coda_window = stream.slice(dt_Coda, dt_Coda+SD)[0].data
    except IndexError as e:
        print("\n %s" % e)
        return
    
    if not (coda_window.size == SD*sr + 1):
        if not (coda_window.size == SD*sr):
            if not (coda_window.size == SD*sr -1):
                if Settings["verbose"]:
                    print("\n\tSize of Coda-window is %s (must be %s)" % (coda_window.size, (SD*sr + 1)))
                # just warning...
                return
    
    #===========================================================================
    # *** Crucial place for calculations ***
    # Add Qp calculation
    
    # RMS of P-bulk-window
    if Pwindow is not None:
        P_window_result = Pwindow.max()#RMS(Pwindow)
    else:
        P_window_result = None
    result.append(P_window_result)
    
    #===
    # S-waves
    
    # RMS of S-bulk-window
    #S_window_result = RMS(Swindow) * 2.5
    # не используем RMS, а максимальную амплитуду в группе S-волн
    S_window_result = Swindow.max()
    result.append(S_window_result)
    
    #mark it!
    ax.plot(seconds_S, S_window_result, "or", markersize=12, markeredgecolor="w", 
        alpha=0.75, zorder=9999)
    
    # for Coda window (original filtered amplitudes)
    # сумма квадратов на кол-во значений (RMS - root mean square)
    #NO! just calc RMS and for S-waves should * 2 (or 2.5 or 3)
    
    #coda_value = RMS(coda_window) # 2 * RMS???
    coda_value = coda_window.max()
    
    if not maxPeak == coda_value:
        print("\n\tWARNING:")
        print("Must be equal %.5f and %.5f" % (maxPeak, coda_value))
    
    # save maximum value in Coda-window
    result.append(coda_value)
    
    # mark coda value
    ax.plot(seconds_Coda1, coda_value, "oc", markersize=12, markeredgecolor="w", alpha=0.75, zorder=9999)
    
    # save RMS!
    #result.append(coda_value)
    # will uase no RMS(coda), but A40 (or A100 or smth.)
    result.append(y100)
    
    # also save koef of slope k
    result.append(_k)
    #===========================================================================
    
    # plot ec. details
    y_pos = original_trace.data.max() / 2
    
    #=== selected windows
    # plot S-window (Swindow)
    S_time = trace.times()[:Swindow.size] + seconds_S
    ax.plot(S_time, Swindow, "-r", lw=2.)
    
    # Coda (coda_window)
    Coda_times = trace.times()[:coda_window.size] + seconds_Coda1
    ax.plot(Coda_times, coda_window, ":r", lw=1.2)
    
    # plot SNR value text
    ax.text(x=seconds_E+0.5, y=y_pos, s="%.1f" % SNR) # SNR
    
    # limits
    ax.set_xlim(seconds_E - 2, seconds_Coda2 + 2)
    ax1.set_ylim(original_trace.data.min(), original_trace.data.max())
    ax.set_ylim(0, trace.data.max()*1.5)
    
    ax.legend()
    
    #plt.tight_layout()
    
    if PLOT:
        plt.show()
    else:
        # just save Figure
        save_to_dir = os.path.join("img", STATION, str(stream[0].stats.starttime)[:19].replace(":", "-"))
        if not os.path.exists(save_to_dir): os.makedirs(save_to_dir)
        outfilename = os.path.join(save_to_dir, "{}_{}_{}__{}s.png".format(STATION, Settings["channel"], Freq, SD))
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
    parser.add_argument("--sd", type=int, 
        help="Specify length of Coda-window (in seconds)")
    
    # maybe we want to plot anyway
    parser.add_argument('--plot', help='Plot results', action='store_true')
    
    # verbose mode
    parser.add_argument("-v", "--verbose", help="output verbosity", 
        action="store_true")
    
    args = parser.parse_args(args)
    
    return args


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
        Settings['channel'] = args.channel[0]#  or Settings['channel']
    
    # may want another SD:
    if args.sd is not None:
        Settings["sd"] = args.sd
        SD = args.sd
    
    # verbose or not
    Settings["verbose"] = True if args.verbose else False
    
    if Settings["verbose"]:
        for k, v in Settings.items():
            print(k, v)
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
    # ======
    
    # get name if input file
    STATION = Settings["station"]
    
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
    
    # 2 * Ts (maybe 2.5*Ts or 1.5*Ts or even 3*Ts)
    KOEF = Settings["koef"] # KOEF
    
    # ABSolute time for start of Coda-window or 2 * Ts or smth
    ABS = False#Settings["absolute"]
    
    # what data to fetch (we will use only 1 channel)
    #which_channel = "?H%s" % (Settings["channel"] if not Settings["rms"] else "?")
    which_channel = "%s" % (Settings["channel"] if not Settings["rms"] else "?")
    
    # NETWORK station etc
    #net, sta, loc, cha = "BR", STATION, "00", which_channel
    # network, station, lcoation, channel
    net, sta, loc, cha = "XQ", STATION, "", which_channel
    
    
    # read data Stream
    DATAFILE = os.path.join("data", STATION, DATAFILE)
    # will cache overall data in memory? read overall data onto memory
    STREAM = read_mseed(DATAFILE, format="MSEED")
    
    # input file in dir "catalog/[STA]/"
    InputFile = os.path.join("catalogs", STATION, InputFile)
    if not os.path.exists(InputFile):
        print("Input file '%s' not found!" % InputFile)
    
    # load input file
    data = np.genfromtxt(InputFile, dtype=None, autostrip=True, names=True,
        loose=True, invalid_raise=False, encoding="utf-8")
    # имена столбцов
    names = data.dtype.names
    # проверить что все нужные столбцы в файле есть
    for col in ("DATE_E", "TIME_E", "TIME_P", "TIME_S"):
        if not col in names:
            print("Coulmn '{}' not found in input file. Fix it!".format(col))
            sys.exit(1)
    # write column headers: DIST_KM AZIMUTH
    names = list(names) + ["DIST_KM", "SNR", "AZIMUTH_AB", "AZIMUTH_BA"]

    #===========================================================================
    # start iterating over Catalog
    
    # save final results
    R = []
    
    for nnn, line in enumerate(data):
        # nice output
        s = "%s of %s\t%s" % (nnn, data.size, " ".join(map(str, line))[:37])
        sys.stdout.write("\r" + s)
        sys.stdout.flush()
        # обрабатывать данные по событию
        parts = dict(zip(names, line))
        
        # check Energy
        ##########if parts["K"] < 9: continue
        
        # считать расстояние, азимуты от станции до события
        lat, lon = parts["LAT"], parts["LON"]
        
        # gps2dist_azimuth(lat1, lon1, lat2, lon2)
        try:
            dist, azAB, azBA = gps2dist_azimuth(Settings["station_lat"],
                Settings["station_lon"], lat, lon)
        except ValueError as e:
            print(e)
            dist = 9999
        # Distance must be in kilometers (not in m)
        dist /= 1000.
        
        # Distance checking...
        #if dist > 500: continue
        
        # where the Coda will start..
        Tcoda = Settings["Tcoda"]
        
        #== make UTC datetime(s)
        
        # time in source
        DT = UTCDateTime(str(parts["DATE_E"]) + u"T" + str(parts["TIME_E"]), precision=2)
        #DT = UTCDateTime(unicode(parts["DATE_E"]) + u"T" + unicode(parts["TIME_E"]), precision=2)
        
        # time-P maybe `None` somewhere... at STDB old Geon data
        # just skip ip
        if parts["TIME_P"].upper() == "NONE":
            continue
        else:
            DT_P = UTCDateTime(parts["DATE_E"] + "T" + parts["TIME_P"], precision=2)
        
        # may be no value for S-obset
        if parts["TIME_S"] == "None":
            # calc from time-P
            DT_S = DT + (DT_P - DT) * CONST_Tp_Ts
        else:
            DT_S = UTCDateTime(parts["DATE_E"] + "T" + parts["TIME_S"], precision=2)

        # check time valid
        assert DT <= DT_P <= DT_S, "Error (%s <= %s <= %s) 'DT <= DT_P <= DT_S' for %s" % (DT, DT_P, DT_S, line)
        
        #=== Fetch data: request this waveforms from MSEED volume (in seconds)
        
        if ABS:
            SEISMOGRAMM_MAX_LEN = Tcoda + SD + 1
        else:
            # calc `Tcoda` as 2*Ts
            SEISMOGRAMM_MAX_LEN = (DT_S - DT) * KOEF + SD + 1 # KOEF == 2 normally
        
        # we will fetch data for this times:
        REQUEST_WINDOW = (0, SEISMOGRAMM_MAX_LEN) # after: Coda starts at 40s + 40 = 120
        
        # Retrieve data
        t1, t2 = DT + REQUEST_WINDOW[0], DT + REQUEST_WINDOW[1]
        
        
        kwargs2 = {"stream": STREAM, 'network': net, 'station': sta, 'location': loc,
            'channel': cha, 'starttime': t1, 'endtime': t2, 'event': 'event'}
        
        # Fetch data
        stream = get_waveforms(**kwargs2)
        
        # maybe duplicate files?
        if stream.count() == 6: stream.merge(method=1)
        
        if stream.count() == 0:
            if Settings["verbose"]: print("\nNo data for event `%s`" % DT)
            continue
        
        # must be all 3 components
        if Settings["rms"] and not stream.count() == 3:
            print("\nERROR: must be 3 channel to make RMS compomnent!")
            continue
        
        
        # SMTH stupid when working with IRIS 20 Hz data TLY
        if stream.count() > 1:
            _starttimes = [tr.stats.starttime for tr in stream]
            #assert _starttimes[0] == _starttimes[1] == _starttimes[2], "stream: %s" % stream
            if not (_starttimes[0] == _starttimes[1] == _starttimes[2]):
                print("stream `%s` is not aligned!" % stream)
                continue
            # also check endtimes!!!
            #_endtimes = [tr.stats.endtime for tr in stream]
            #assert _endtimes[0] == _endtimes[1] == _endtimes[2]
        
        # detrend (this includes removing mead value)
        stream.detrend('demean')
        
        #=== Remove response first, then make RMS componet etc
        # remove response or not
        if Settings["simulate"]:
            try:
                stream = remove_response(stream, STATION)
            except ValueError as e:
                print("\tERROR: remove response failed: %s. Sream `%s`" % (e, stream))
                raw_input()
                continue
        
        #===
        
        if Settings["rms"]:
            # **2
            __N = np.power(stream.select(channel='N')[0].data, 2)
            __E = np.power(stream.select(channel='E')[0].data, 2)
            __Z = np.power(stream.select(channel='Z')[0].data, 2)
            # R - result of calculations
            RMS_component = np.sqrt( (__N + __E + __Z) / 3 )
            # remove 1st 2 components
            for _ in range(2): stream.pop(0)
            stream[0].data = RMS_component
            stream[0].stats.channel = "RMS"
            #===
            
        assert stream.count() == 1, "Must be 1 trace! we have: %s" % stream
        
        #=============
        #== считать окна, для времен вступлений, 1 раз - для всех частот
        # начало и конец файла
        T0, Tend = stream[0].stats.starttime, stream[0].stats.endtime
        
        # this check is useless:
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
        
        # ! check we will make: T0 must be EXACTLY what we fetched...
        if not T0 == t1:
            # first check the difference...
            diff = abs(T0 - t1)
            # if difference < 1 second
            if diff < 1: pass
            else:
                if Settings["verbose"]: 
                    print("\nStream must start at the same time as T1 (%s != %s)!!! \nDiff = %s" % (T0, t1, diff))
                continue
        
        #=== Время в очаге T0: сколько это в секундах относительно начала файла?
        # время в очаге
        seconds_E = DT - T0
        # P-wave onset
        seconds_P = DT_P - T0
        # S-wave onset
        seconds_S = DT_S - T0
        
        #===
        # WARNING! recalc distance according to formula 8 * (Ts - Tp + 0.2)
        ##real_dist = 8 * (seconds_S - seconds_P + 0.2)
        ##if Settings["verbose"]:
        ##    print("\nReal_dist = %.2f, dist was %.2f" % (real_dist, dist))
        ##dist = real_dist
        
        #===============================================================
        # для всех частот считать...
        # порядок записи частот: частота порог
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
                    Tcoda
                )
                plt.close()
            #except BaseException as e:
            #    print("!"*33)
            #    print(e)
            #    print("!"*33)
            #    print()
            #    continue
            #    
            except KeyboardInterrupt:
                print("Interrupted by user")
                PLOT = False
                continue
            
            #===========================================================
            # if no values return, skip this freq...
            if (result is None) or ( len(result) == 1 ):
                if Settings["plot"]: plt.close()
                # надо считать по другим частотам, если по 1-й неудачно? Yeas!
                calc_freq_values.append( 0 )
                continue
            
            # parse results
            # SNR, P-value, S-value, Coda-value
            SNR, P_value, S_value, C_value, Y100, _k = result

            #====================
            # считать логарифм отношения Амплитуды прямой волны (с коррекцией на геом. расхождение) и коды 
            
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # most important calculation: LN(As * Z(R) / Ac)
            #########LN_As_Ac = np.log( P_Z / (C_Z * Z_R ) )
            #LN_Ap_Ac = np.log( P_value / C_value )
            LN_Ap_Ac = np.log( P_value / Y100 )
            
            #LN_As_Ac = np.log( S_value / C_value )
            LN_As_Ac = np.log( S_value / Y100 )
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            calc_freq_values.append( LN_As_Ac )
            
            # save distance-Amplitude - Ap/As/Ac/Y100 values
            datapath = os.path.join("values", STATION)
            if not os.path.exists(datapath): os.makedirs(datapath)
            # filename: STA_FREQ_CH_vals.txt
            datafilename = "{0}_{1:0.2f}_{2}_vals_SD{2}__ABS_{3}_T{4}.txt".format(STATION, Freq, SD, ABS, 0 if Tcoda is None else Tcoda)
            datafilename = os.path.join(datapath, datafilename)
            if not os.path.exists(datafilename):
                # create datafile and add Header(s)
                _dataf = open(datafilename, "w")
                _headers = "\t".join("DATE_E TIME_E LAT LON K CHA DIST AZab SNR Ap As Ac Y100 LN(Ap/Ac100) LN(As/Ac100) _k".split())
                _dataf.write(_headers + "\n")
                _dataf.close()
            
            with open(datafilename, "a") as _dataf:
                # write idE, idDir, DATE_E, TIME_E, LAT, LON, K, CHA, SNR, P_value, S_value_ C_value, C_Z, LN
                s = "{DATE_E}\t{TIME_E}\t{LAT:.3f}\t{LON:.3f}\t{K:.1f}\t".format(**parts)
                s += "{0}\t".format(Settings["channel"])# Channel Name
                # distance: real and for calc. + Azimuth
                s += "%.1f\t%.1f\t" % (dist, azAB)
                # append results (SNR, P_value, S_value, C_value, A100)
                s += "\t".join(["%.5f"%_r for _r in result])
                
                # append final resultL LN(As/Ac)
                s += "\t%.4f\t%.4f" % (LN_Ap_Ac, LN_As_Ac)
                
                # slope of coda exponent _k
                s += "\t%.4f" % _k
                
                _dataf.write(s + "\n")
                _dataf.close()

            #=== calculations ended
            #===========================================================
            
        # сохранить пары значений (расстояние_в_км, Логарифм_отношения_по_P)
        R.append( {dist: calc_freq_values} )
    
    # finished for all events from Catalog
    #===========================================================================
    
    print()
    print("Results:")
    print(":"*33)
    
    # now lets calc further: для каждой частоты, строить график
    for freq_num, Freq in enumerate(Freqs):
        # remove null values in list
        R = [r for r in R if r.values() != [[]] ]
        print()
        print('R')
        print(R)
        print()
        print()
        
        # get values for x and y
        
        #x = np.array( [r.keys()[0] for r in R] )
        #y = np.array( [r.values()[0][freq_num] for r in R] )
        
        x, y = [], []
        for r in R:
            for k, v in r.items():
                x += [k]
                y += [v[freq_num]]
        
        # кое-где у нас 0 в y - убрать соотв x
        _indices = np.where(y != 0)
        y = y[_indices]
        x = x[_indices]
        
        # calc coeff correlation (after removing zeros 0)
        koef_pears, p = pearsonr(x, y)
        
        # calc regression (polyfit)
        try:
            m, b = np.polyfit(x, y, 1)
        except TypeError:
            print("linear fitting error: %s (frequency=%.2f Hz). Skipping..." % (e, Freq))
            continue
        Y = m * x + b
        
        # calc by ROBUST method
        try:
            a, b2 = linear_fit(y, x)
        except BaseException as e:
            print("robust fitting error: %s" % e)
            a, b2 = 1, 0
        Y2 = a * x + b2
        
        # plot Linear regr. x-y
        fig = plt.figure()#dpi=100, figsize=(16, 9))
        ax = fig.add_subplot(111)
        ax.set_title("Freq = {}".format(Freq))
        ax.plot(x, y, 'ro', label="x-y (Corr=%.3f)" % koef_pears, markeredgecolor="k", markeredgewidth=0.5)
        ax.plot(x, Y, "-k", label=r"$y=%.5f x + %.5f$" % (m, b))
        
        # plot ROBUST results
        ax.plot(x, Y2, "-b", lw=2.,
            label=r"$y=%.5f x + %.5f$" % (a, b2))
        
        ax.legend()
        
        # Далее рассчитать сами значения добротности Q
        V = Settings['vs']# if CALCQ == "P" else Settings['vs']
        assert V == 3.51
        Q = -1 * np.pi * Freq / (V * m)
        
        # 2nd method - Robust mean
        Q2 = -1 * np.pi * Freq / (V * a)
        
        # insteadn of printing, save these values
        outpfname = "Qvalues_{}_SD{}.txt".format(STATION, SD)
        if not os.path.exists(outpfname):
            with open(outpfname, "w") as _f:
                # header:
                # FREQ  Q RobustQ  Corr N
                _f.write('\t'.join("CHANNEL FREQ Q RobustQ Corr N".split()) + "\n")

        _f = open(outpfname, "a")
        
        # write values
        s = '%s\t%.2f\t%.1f\t%.1f\t%.3f\t%d\n' % (Settings["channel"], Freq, Q, Q2, koef_pears, x.size)
        _f.write(s)
        print(s.strip())
        _f.close()
        
        # save to img/STA/ folder
        save_to_dir = os.path.join("img", STATION)
        if not os.path.exists(save_to_dir): os.makedirs(save_to_dir)
        outfilename = os.path.join(save_to_dir, "{0:02d}__{1}_{2}_{3}__{4}s__ABS-{5}_T{6}.png".format(freq_num+1, STATION, Settings["channel"], int(Freq), SD, ABS, 0 if Tcoda is None else Tcoda))
        plt.savefig(outfilename)
        #plt.show()
        plt.close()

