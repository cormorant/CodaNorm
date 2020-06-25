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
Программа расчета среднеквадратических амплитуд прямых волн и коды сейсмограмм.
Входной формат данных - GSE2 (опционально, добавить MSEED, исходный Байкал-5).

Расчет осуществляется по 2-м алгоритмам:
    - алгоритм нормализации.
    - по амплитудным спектрам коды, волн P, S.
UPDATE:
    - добавлен(?) метод огибающей коды

UPDATE 2020-06-25:
    - added command line arguments, in order to переопределять "coda.conf" params
    
"""
APP_NAME = "CodaNorm"
__version__="0.2.0.1"
COMPANY_NAME = 'GIN SB RAS'

import os
import sys
import time
import datetime
#import math

import obspy
from obspy import UTCDateTime, Trace, Stream, read as read_mseed, read_inventory
from obspy.geodetics.base import gps2dist_azimuth
#from obspy.signal.invsim import simulate_seismometer, corn_freq_2_paz
#for evelope calculations
import obspy.signal

import numpy as np
from scipy.stats.stats import pearsonr

# CodaNormlib code
from CodaNormLib import (STATIONS, linear_fit, read_settings, RMS, 
    calc_signal_noise_ratio, get_waveforms)

#import json

# command line arguments may override default settings from `coda.conf` file
import argparse

# grapihcs
import matplotlib.pyplot as plt
#from matplotlib import mlab

#from itertools import cycle
#COLORS = cycle(("g", "c", "r"))

# This doesnt work:
#import warnings
#def fxn():
#    warnings.warn("deprecated", FutureWarning)#
#
#with warnings.catch_warnings():
#    warnings.simplefilter("ignore")
#    fxn()


CHANNELS = ("N", "E", "Z")

CONST_Tp_Ts = 1.7


def module_path():
    ''' add link to where this coda is stolen '''
    if hasattr(sys, "frozen"):
        return os.path.dirname(sys.executable)
        #return os.path.dirname(unicode(sys.executable, sys.getfilesystemencoding( )))
    return os.path.dirname(__file__)
    #return os.path.dirname(unicode(__file__, sys.getfilesystemencoding( )))


# get current dir, may vary if run from EXE file
CurrDir = module_path()

CONFIG_FILENAME = os.path.join(CurrDir, "coda.conf")

#===============================================================================


def get_Tabs_and_alpha(dist):
    # we will Hard-code Tabs:
    # R           2Ts max = 
    #   0- 50 km: 33
    #  50-100 km: 56 
    # 100-150 km: 80
    # 150-200 km: 104
    # 200-250 km: 128
    
    if dist <= 50:
        Tcoda = 33
    elif 50 < dist <= 100:
        Tcoda = 56
    elif 100 < dist <= 150:
        Tcoda = 80
    elif 150 < dist <= 200:
        Tcoda = 104
    elif 200 < dist <= 250:
        Tcoda = 128
    else:
        raise NotImplementedError("Must choose DIST to calc from range '1-3'")
     
    #! also Geom. spreadin will change:
    #ALPHA = Settings["alpha"]
    # Alpha - parameter of geom. spreading, depends on Distance we calc:
    # for R <= 100, alpha = -1 (but for R between 60 and 100, R always == 60)
    # for R > 100, alpha = -0.5
    
    #ALPHA = -0.5 if dist == 3 else -1
    
    # alpha = -0.5 also for low freqs, and for all with dist > 100km
    #ALPHA = -0.5 if dist > 100 else -1
    ALPHA = ALPHA2 if dist > 100 else ALPHA1
    
    return Tcoda, ALPHA



def calculate(Freq, f1, f2, stream, 
    seconds_E, seconds_P, seconds_S, dt, 
    dt_P, dt_S, Settings, dist, Tcoda):#azBA=0):
    """ 
    Freq, f1, f2
    stream - a copy of 'slices' original STREAM object that saves all data
    seconds_E - сколько это в секундах относительно начала файла - время в очаге
    (всегда 0 если у нас REQUEST_WINDOW = (0, ...))
    
    seconds_S - сколько прошло до S-волны в секундах от начала файла (или T0)
    
    расчет параметров коды - по горизонтальным компонентам """
    
    # save original trace
    original_trace = stream[0].copy()
    T0 = original_trace.stats.starttime
    
    #sr = original_trace.stats.sampling_rate # sampling rate, 100 Hz normally
    sr = float("%.1f" % original_trace.stats.sampling_rate) #Why 1 digit after .
    
    # stream must start at T0
    if abs(seconds_E - 0.) > 1:
        assert seconds_E == 0. , "\nseconds_E '%s' !=5" % (seconds_E)
    # time of Coda
    
    # according to [Aki 1980] [Pavlenko 2008] [Indus etc.]
    # нужно брать коду за фиксированныое время - например 40 с после Т0 (время в очаге)
    # однако для расстояний > 70 км, 40 сек будет раньше времени 2*Ts !!!
    # поэтому по формуле нужно пересчитать амплитуду коды (например на 80 с)
    # в амплитуду коды на 40с
    # по формуле
    # RMS(Coda at 80s) * D(f, 40) / D(f, 80)
    # где D(f, t) - интенсивность кода-волн со временем t
    
    # real 2*Ts
    dt_Coda_shouldbe = dt + 2. * (dt_S - dt)
    
    if not ABS:
        # 2 * tS
        dt_Coda = dt + Settings["koef"] * (dt_S - dt)
    # время в очаге + 50 с ???
    else:
        # `dt_Coda` - is UTCDatetime obj that saves Absolute time of Coda-window
        dt_Coda = dt + Tcoda
    
    #dt_Coda = dt_Coda_shouldbe
    # seconds
    seconds_Coda1 = dt_Coda - T0
    seconds_Coda2 = seconds_Coda1 + SD
    
    # IMPORTANT UPDATE:
    # calc SD1 manually (if not set)
    SD1 = Settings['sd1']
    if SD1 == 0:
        # window for Direct body-waves = 0.25 * Ts
        SD1 = int( round(0.25*seconds_S) )
        
        SD1P = int(round(SD1 / 2))
    
    # check: S-window and Coda-window may intersects!
    if (seconds_Coda1 - seconds_S) < SD1:
        print("\nWindow S and Coda intersects for %s!!!" % dt)
        # save data aboutIntersections
        '''
        _f = open("intersects.txt", "a")
        s = "%s\t%s (P=%s, S=%s)\n" % (STATION, dt, dt_P, dt_S)
        _f.write(s)
        _f.close()
        '''
        return

    # get results for this filename at this freqs
    if PLOT:
        fig, ax = plt.subplots(figsize=(16, 9), nrows=1, sharex=True)
        fig.suptitle("Dist = {:.0f} \tStream {}--{}. Freq {} Hz ({}-{})".format(dist, 
            stream[0].stats.starttime, stream[0].stats.endtime, Freq, f1, f2), fontsize=14)
    else:
        axes = np.arange(stream.count())
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
    
    
    if trace.data.max() == trace.data.min():
        print("\nWarning: may we have here solid lines!!!!!!!!")
        return
    
    # check SNR on filtered signal!
    # noise windows will always be = 3 second
    # calc at the end of the Coda-window!
    SNR = calc_signal_noise_ratio(stream, dt, dt_Coda+SD, 3)
    #SNR = calc_signal_noise_ratio(stream, dt, dt_Coda, 5)
    
    result += [SNR]
    
    if SNR < Settings['minsnr']:
        if Settings["verbose"]: 
            print("\nSNR for freq %s is too SMALL (%.1f)!..." % (Freq, SNR))
        return result
    
    if PLOT:
        color = 'g'
        # original signal
        ax.plot(original_trace.times(), original_trace.data, color, 
            label=trace.stats.channel, lw=.25, zorder=111)
        # plot filtered
        ax.plot(trace.times(), trace.data, "k", lw=.5, zorder=222)
        
        # etc
        
        # mark time of Event, P and S time by vertical lines, start -- end
        ax.axvline(x=seconds_E, linestyle="--", color="y", lw=2.) # Event
        ax.axvline(x=seconds_P, linestyle="--", color="m", lw=0.25) # P
        ax.axvline(x=seconds_P+SD1P, linestyle="--", color="m", lw=0.25) # P+SD1p
        
        # S-window
        ax.axvline(x=seconds_S, linestyle="--", color="k", lw=1) # S
        ax.axvline(x=seconds_S+SD1, linestyle="--", color="k", lw=1) # S+SD1
        # mark coda
        ax.axvline(x=seconds_Coda1, linestyle="-", color="r", lw=2.) # coda
        ax.axvline(x=seconds_Coda2, linestyle="-", color="r", lw=2) # coda
        
        # anyway plot 2*TS marker
        seconds_Coda_2ts = dt_Coda_shouldbe - T0
        ax.axvline(x=seconds_Coda_2ts, linestyle=":", color="c", lw=2.5) # coda 2TS
    
    #===================================================================
    #=== все вычисления по данной компоненте
    
    # получить сами массивы данных (отфильтрованных)
    
    # S-window data
    Swindow = stream.slice(dt_S, dt_S+SD1)[0].data
    
    # Coda-window data
    coda_window = stream.slice(dt_Coda, dt_Coda+SD)[0].data
    
    if not (coda_window.size == SD*sr + 1):
        if not (coda_window.size == SD*sr):
            if not (coda_window.size == SD*sr -1):
                print("Size of Coda-window is %s (must be = %s)" % (coda_window.size, (SD*sr + 1)))
                # just warning...
                #return
    
    #=== Calculations (normally on the Envelope)
    
    # uses "ENVELOPE"
    
    # envelope of S-bulk window (Hilbert transform)
    valuesS = obspy.signal.filter.envelope(Swindow)
    
    # original values for Coda
    valuesC = coda_window
    
    #===========================================================================
    # *** Crucial place for calculations ***
    
    # just Max-Min (may lead to Very bad results and too low correlation!!!)
    S_window_result = valuesS.max() - valuesS.min()
    result.append(S_window_result)
    
    # for Coda window (original filtered amplitudes)
    # сумма квадратов на кол-во значений (RMS - root mean square)
    #NO! just calc RMS and for S-waves should * 2 (or 2.5 or 3)
    coda_value = RMS(valuesC) # 2 * RMS
    #coda_value = valuesC.max()
    
    result.append(coda_value)
    #===========================================================================
    
    # plot details
    if PLOT:
        y_pos = original_trace.data.max() / 2
        
        #=== envelopes
        # plot envelope of S-window (Swindow)
        S_time = trace.times()[:Swindow.size] + seconds_S
        ax.plot(S_time, valuesS, "--r", lw=2.)
        
        # Coda envelope (coda_window)
        #Coda_times = trace.times()[:coda_window.size] + ABS
        Coda_times = trace.times()[:coda_window.size] + seconds_Coda1
        ax.plot(Coda_times, coda_window, "--r", lw=2.)
        
        #=== calc values
        # calc and mark MAX-MIN of envelope
        max_min_bulk = Swindow.max() - Swindow.min()
        
        # coda amplitude:
        coda_amplitude = coda_window.max() - coda_window.min()
        
        #=== Texts:
        # plot SNR value
        ax.text(x=seconds_E+2, y=y_pos, s="%.1f" % SNR) # SNR
        
        # plot RMS value of bulk-window
        ax.text(x=seconds_S+SD1/2, y=y_pos, s="RMS(bulk)=%.1f (Max-min=%.1f)" % (result[1], max_min_bulk)) # RMS(Swindow)
        
        # plot RMS value of Coda-window
        ax.text(x=seconds_Coda2-10, y=y_pos, s="RMS(coda)=%.1f (Max-min=%.1f)" % (result[2], coda_amplitude)) # RMS(Coda_window)
        
        
        ax.set_xlim(seconds_E - 5, seconds_Coda2 + 5)
        ax.set_ylim(original_trace.data.min(), original_trace.data.max())
        
        ax.legend()

        plt.show()
        '''
        save_to_dir = os.path.join("img", str(stream[0].stats.starttime)[:19].replace(":", "-"))
        if not os.path.exists(save_to_dir): os.makedirs(save_to_dir)
        outfilename = os.path.join(save_to_dir, "{}_{}.png".format(STATION, Freq))
        plt.savefig(outfilename)
        '''
        plt.close()
    # the end
    return result


def parse_command_line_args(args=None):
    parser = argparse.ArgumentParser(description='CodaNorm arguments that override default settings')
    
    # argument on which channel to calc
    parser.add_argument("-c", "--channel", help="Channel to use (N, E, Z).", 
        choices=CHANNELS, nargs='+')
    
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
    
    # verbose or not
    Settings["verbose"] = True if args.verbose else False
    
    if Settings["verbose"]:
        for k, v in Settings.items():
            print(k, v)
        # nice wait
        _N = 10
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
    # ======
    
    # get name if input file
    STATION = Settings["station"]
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
    ABS = Settings["absolute"]
    
    # geomettrical spreading params:
    R1 = Settings["r1"]
    R2 = Settings["r2"]
    ALPHA1 = Settings["alpha1"]
    ALPHA2 = Settings["alpha2"]
    
    # what data to fetch (we will use only 1 channel)
    net, sta, loc, cha = "BR", STATION, "00", "?H%s" % Settings["channel"]
    
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
        if dist > 250: continue
        
        #TODO: no need to hard-code 60 and 100; use vars R1 and R2 instead
        real_dist = 0# in order not to use same VAR as dist
        real_dist = dist

        # if less then 60 km
        if R1 < dist <= R2:
            # если расстояния от 60 до 100 км, геом расхождение будет 60 ^ -1
            dist = R1
        
        # get Tabs and ALPHA param for this distance:
        Tcoda, ALPHA = get_Tabs_and_alpha(real_dist)
        
        # make UTC datetime(s)
        DT = UTCDateTime(parts["DATE_E"] + "T" + parts["TIME_E"], precision=2)
        DT_P = UTCDateTime(parts["DATE_E"] + "T" + parts["TIME_P"], precision=2)
        # may be no value for S-obset
        if parts["TIME_S"] == "None":
            # calc from time-P
            DT_S = DT + (DT_P - DT) * CONST_Tp_Ts
        else:
            DT_S = UTCDateTime(parts["DATE_E"] + "T" + parts["TIME_S"], precision=2)
        

        #WARNING: for UZR (A-7b) move P-S-onsets window -0.5s
        DT_P = DT_P - 0.25
        DT_S = DT_S - 0.25
        # check time valid
        assert DT <= DT_P <= DT_S, "Error (%s <= %s <= %s) 'DT <= DT_P <= DT_S' for %s" % (DT, DT_P, DT_S, line)
        
        #=== Fetch data
        # request this data from MSEED volume (in seconds)
        # full length of seismogramm no more then Tabs + Coda_length (for ex.: 40s+40s = 80s)
        SEISMOGRAMM_MAX_LEN = Tcoda + SD + 1
        REQUEST_WINDOW = (0, SEISMOGRAMM_MAX_LEN) # after: Coda starts at 40s + 40 = 120
        
        # Retrieve data
        t1, t2 = DT + REQUEST_WINDOW[0], DT + REQUEST_WINDOW[1]
        
        
        kwargs2 = {"stream": STREAM, 'network': net, 'station': sta, 'location': loc,
            'channel': cha, 'starttime': t1, 'endtime': t2, 'event': 'event'}
        
        # Fetch data
        stream = get_waveforms(**kwargs2)
        
        if stream.count() == 0:
            if Settings["verbose"]: print("\nData missing for event %s" % DT)
            #_f = open(r"d:\Work\python\_new\CodaNorm\Missing.txt", "a")
            #_f.write("%s\n" % DT)
            #_f.close()
            continue
        
        # detrend (this includes removing mead value)
        stream.detrend("linear")
        
        # simulate or not:
        if Settings["simulate"]:
            # for TLY
            if STATION == "TLY":
                inv = read_inventory("resp/IRIS__TLY__2003_12.xml")
                stream.remove_response(inventory=inv, output="VEL", plot=False)
            elif STATION == "MXMB":
                inv = read_inventory("resp/BR__MXMB_7hr.xml")
                # pre-filter outside of desired frequency range
                pre_filt = [0.01, 0.1, 30, 50]
                stream.remove_response(inventory=inv, output="VEL", plot=False,
                    pre_filt=pre_filt)
            elif STATION == "UUDB":
                respfile = "resp/BR_UUDB__Centaur.resp"
                inv = read_inventory(respfile)
                stream.remove_response(inventory=inv, output="VEL", plot=False)
            elif STATION == "KELR":
                respfile = "resp/BR_KELR__CMG-40T.resp"
                inv = read_inventory(respfile)
                pre_filt = [0.01, 0.1, 20, 25]
                stream.remove_response(inventory=inv, output="VEL", plot=False, pre_filt=pre_filt)
            else: pass
            # do it
            #stream.remove_response(inventory=inv, output="VEL", plot=False)

        #== считать окна, для времен вступлений, 1 раз - для всех частот
        # начало и конец файла
        T0, Tend = stream[0].stats.starttime, stream[0].stats.endtime
         
        if not T0 == t1:
            # first check the difference...
            diff = abs(T0 - t1)
            # if difference < 1 second
            if diff < 1: pass
            else:
                if Settings["verbose"]: 
                    print("\nStream must start at the same time as T1 (%s != %s)!!! \nDiff = %s" % (T0, t1, diff))
                continue
        
        # also check EndTime
        if not Tend == t2:
            # check the difference...
            diff = abs(Tend - t2)
            if diff < 1: pass
            else:
                print("\nStream must end at the same time as T2 (%s != %s)!!! \nDiff = %s" % (Tend, t2, diff))
                continue
        
        # check all traces start at the same time
        for tr in stream:
            # check diff
            diff = T0 - tr.stats.starttime
            if diff > 0:
                print("T0 != tr.stats.starttime : (%s and %s)" % (T0, tr.stats.starttime))
                print("Diff = %s" % diff)
                raw_input()
                continue
        
        # но могут быть и 2 сейсмограммы. выбрать более длинную или 1st?
        if stream.count() > 1:
            #print("\nUsing FIRST stream instead!!!")
            stream = stream[:1]
        
        #=== Время в очаге T0: сколько это в секундах относительно начала файла?
        # время в очаге
        seconds_E = DT - T0
        # P-wave onset
        seconds_P = DT_P - T0
        # S-wave onset
        seconds_S = DT_S - T0
        
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
                    Settings, real_dist or dist,#if real_dist == 0 will use dist
                    Tcoda
                )
            #except BaseException as e:
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
            SNR, P_Z, C_Z = result

            #====================
            # считать логарифм отношения Амплитуды прямой волны (с коррекцией на геом. расхождение) и коды 
            
            # also for distance R <60 and 60-100km, we use ALPHA = -1: R^-1 (and 60^-1 for 60-100km)
            
            # geom. spreading Z(R)
            Z_R = np.power(float(dist), ALPHA)
            
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # most important calculation: LN(As * Z(R) / Ac)
            LN_As_Ac = np.log( P_Z / (C_Z * Z_R ) )
            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
            calc_freq_values.append( LN_As_Ac )
            
            # save distance-Amplitude (LN(As/Ac with ALPHA correction) values
            datapath = os.path.join("values", STATION)
            if not os.path.exists(datapath): os.makedirs(datapath)
            # filename: STA_FREQ_CH_vals.txt
            datafilename = "{0}_{1:0.2f}_{2}_vals_SD{3}.txt".format(STATION, Freq, Settings["channel"], SD)
            datafilename = os.path.join(datapath, datafilename)
            if not os.path.exists(datafilename):
                # create datafile and add Header(s)
                _dataf = open(datafilename, "w")
                _headers = "\t".join("DATE_E TIME_E LAT LON K CHA ALPHA DIST REAL_DIST SNR P_Z C_Z LN".split())
                _dataf.write(_headers + "\n")
                _dataf.close()
            

            with open(datafilename, "a") as _dataf:
                # write idE, idDir, DATE_E, TIME_E, LAT, LON, K, CHA, ALPHA, SNR, P_Z, C_Z, LN
                s = "{DATE_E}\t{TIME_E}\t{LAT:.3f}\t{LON:.3f}\t{K:.1f}\t".format(**parts)
                s += "{0}\t{1}\t".format(Settings["channel"], ALPHA)
                # distance: real and for calc.
                s += "%.1f\t%.1f\t" % (dist, real_dist)
                # append results (SNR, P_Z, C_Z)
                s += "\t".join(["%.3f"%_r for _r in result])
                # append final resultL LN(As/Ac with Alpha correction)
                s += "\t%.3f" % LN_As_Ac
                
                _dataf.write(s + "\n")
                _dataf.close()

            #=== calculations ended
            #===========================================================
            
        # after calculations, restore real distance value
        dist = real_dist
        
        # сохранить пары значений (расстояние_в_км, Логарифм_отношения_по_P)
        R.append( {dist: calc_freq_values} )
    
    # finished for all events from Catalog
    #===========================================================================
    
    # stop for now
    sys.exit(0)
    
    print()
    print("Results:")
    print(":"*33)
    
    # now lets calc further: для каждой частоты, строить график
    for freq_num, Freq in enumerate(Freqs):
        # remove null values in list
        R = [r for r in R if r.values() != [[]] ]
        # get values
        x = np.array( [r.keys()[0] for r in R] )
        # for freq #1
        y = np.array( [r.values()[0][freq_num] for r in R] )
        
        # кое-где у нас 0 в y - убрать соотв x
        _indices = np.where(y != 0)
        y = y[_indices]
        x = x[_indices]
        
        # calc coeff correlation (after removing zeros 0)
        koef_pears, p = pearsonr(x, y)
        
        # calc regression (polyfit)
        m, b = np.polyfit(x, y, 1)
        Y = m * x + b
        #print x.size,
        
        # calc by ROBUST method
        a, b2 = linear_fit(y, x)
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
        #ax.set_ylim(-5, 0)
        #ax.set_xlim(0, 80)
        
        # Далее рассчитать сами значения добротности Q
        V = Settings['vs']# if CALCQ == "P" else Settings['vs']
        assert V == 3.51
        Q = -1 * np.pi * Freq / (V * m)
        
        # 2nd method - Robust mean
        Q2 = -1 * np.pi * Freq / (V * a)
        
        
        try:
            print('%.1f\t%.1f RobustQ=\t%.1f\tCorr= %.3f\tN= %d' % (Freq, Q, Q2, koef_pears, x.size))
        except TypeError:
            print("'%s' hz freq:\terror calc on..." % Freq)
        
        # save to img/STA/ folder
        save_to_dir = os.path.join("img", STATION)
        if not os.path.exists(save_to_dir): os.makedirs(save_to_dir)
        outfilename = os.path.join(save_to_dir, "{0:02d}__{1}_{2}_{3}.png".format(freq_num+1, STATION, Settings["channel"], int(Freq)))
        #plt.savefig(outfilename)
        plt.show()
        plt.close()

