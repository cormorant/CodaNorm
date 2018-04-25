#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division
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
Body-wave attenuation calculation by coda-normalization method.
Input data format GSE2.
"""
APP_NAME = "CodaNorm"
__version__="0.2"
COMPANY_NAME = 'GIN SB RAS'

import os
import sys
#from obspy.gse2.core import isGSE2, readGSE2
from obspy import read as readGSE2
def isGSE2(): return True

from obspy.core.util.geodetics import gps2DistAzimuth, locations2degrees
from obspy.signal.rotate import rotate_NE_RT

import numpy as np
from scipy.signal import butter, lfilter

import ConfigParser

from itertools import cycle

# lest import lib for writing to excel file
try:
    import xlwt
except ImportError, e:
    print e
    print("Saving to Excel isn't available! Install xlwt library!")
    sys.exit(1)

try:
    import matplotlib.pyplot as plt
except ImportError as e:
    print e
    matplotlib_is_imported = False
else:
    matplotlib_is_imported = True
    COLORS = cycle(("g", "c", "r"))


def calc_spectrum(y, dt):
    """ spectrum calculation  """
    _len = y.size
    half = int(_len/2)#np.ceil
    freqs = np.fft.fftfreq(_len, d=dt)
    freqs = freqs[:half]
    pz = np.fft.fft(y)
    pz = np.fft.fftshift(pz)
    ampsp = np.abs(pz)
    ampsp = ampsp / _len
    if _len % 2 == 0: ampsp = ampsp[half:]
    else: ampsp = ampsp[half:-1]
    ampsp = 2 * ampsp
    return freqs, ampsp


def butter_bandpass_filter(data, lowcut, highcut, fs, order=4):
    ''' bandpass Butterworth filter '''
    nyquist = 0.5 * fs
    low = lowcut / nyquist
    high = highcut / nyquist
    b, a = butter(order, [low, high], btype='band')
    return lfilter(b, a, data)


def calc_seconds_from_time(_timeStr):
    """ calculates seconds since the beginning of the day from a text time """
    hour, minute, second = map(float, _timeStr.split(":"))
    return hour * 3600 + minute * 60 + second


def calc_power(P, S, C):
    """ RMS spectrum estimation.
    Input - spectrum in windows (filtered values P, S, coda).
    """
    return [np.sqrt( np.sum(np.power(a, 2)) / a.size ) for a in (P, S, C)]


def main(Freq, f1, f2, filename, secondsE, secondsP, secondsS, Settings, rotate=None):
    """ calculation itself """
    # plot or not
    PLOT = Settings["plot"]
    if not matplotlib_is_imported: PLOT = False
    # LENGTH OF WINDOWS and KOEF
    SD = Settings["sd"]
    KOEF = Settings["koef"]
    # open ofiginal GSE2 file
    if isGSE2(filename):
        stream = readGSE2(filename)
    else:
        print("File {} is not in GSE2 format!"),
        return
    #===
    result = []
    # get results for this filename at this freqs
    if PLOT:
        fig = plt.figure(dpi=100, figsize=(16, 9))
        ax1 = fig.add_subplot(311)
        ax2 = fig.add_subplot(312, sharex=ax1)
        ax3 = fig.add_subplot(313, sharex=ax1)
        #ax4 = fig.add_subplot(414)
        ax1.set_title("File {}. Freq {} ({}-{})".format(os.path.split(filename)[-1],
            Freq, f1, f2))
        axes = (ax1, ax2, ax3)
    else:
        axes = np.arange(3)
    #=== calculating
    # rotate here traces 1 and 2 (N, E)
    # calculate the spectra of the window P, S, coda
    for ax, trace in zip(axes, stream.traces):
        sr = trace.stats.sampling_rate # sampling rate
        y = butter_bandpass_filter(trace.data, f1, f2, sr)
        # delete mean value
        y -= y.mean()
        if PLOT:
            color = COLORS.next()
            ax.plot(trace.times(), trace.data, color, label=trace.stats.channel, alpha=0.5)
            # plot filtered
            ax.plot(trace.times(), y, "k")
            ax.set_ylim(y.min(), y.max())
            ax.legend()
        #===
        # 0 -- start time of the file
        secondsStart = calc_seconds_from_time(trace.stats.starttime.time.strftime("%H:%M:%S.%f"))
        # if the file starts after the event
        if secondsStart > secondsE:
            print("File {} starts after Event!"),
            return
        # reduce the time in seconds from the beginning of the file
        timeE = secondsE - secondsStart
        timeP = secondsP - secondsStart
        timeS = secondsS - secondsStart
        # calc time of start Coda
        timeCoda = int(np.ceil(KOEF *  (timeS - timeE)) + timeE)
        if PLOT:
            # mark time of Event, P and S time by vertical lines, start -- end
            ax.axvline(x=timeE, linestyle="--", color="y") # Event
            ax.axvline(x=timeP, linestyle="--", color="r") # P
            ax.axvline(x=timeP+SD, linestyle="--", color="r") # P
            ax.axvline(x=timeS, linestyle="--", color="r") # S
            ax.axvline(x=timeS+SD, linestyle="--", color="r") # S
            # mark coda
            ax.axvline(x=timeCoda, linestyle="-", color="r") # coda
            ax.axvline(x=timeCoda+SD, linestyle="-", color="r") # coda
        #=== calculations for component
        #= window for P
        Pindex1 = int( np.ceil(timeP * sr) )
        Pindex2 = int( Pindex1 + SD * sr )
        Pwindow = y[Pindex1:Pindex2] # P
        #= window for S
        Sindex1 = int( np.ceil(timeS*sr) )
        Sindex2 = int( Sindex1 + SD * sr )
        Swindow = y[Sindex1:Sindex2] # S
        # C - кода
        coda_index1 = int(np.ceil(KOEF * (timeS - timeE)) * sr + np.ceil(timeE * sr))
        coda_index2 = int(coda_index1 + SD * sr)
        # coda array
        coda_window = y[coda_index1:coda_index2]
        if not coda_window.size:
            print("Not enough time for coda to cut!"),
            return
        #=== spectrums or just normalize
        # calc RMS values of spectrum
        if "SPECTR" in Settings["algorithm"].upper():
            freqs, valuesP = calc_spectrum(Pwindow, trace.stats.delta) # P
            freqs, valuesS = calc_spectrum(Swindow, trace.stats.delta) # S
            freqs, valuesC = calc_spectrum(coda_window, trace.stats.delta) # C
        # or by amplitudes
        else:
            valuesP, valuesS, valuesC = Pwindow, Swindow, coda_window
        # calc result on windows
        result += calc_power(valuesP, valuesS, valuesC)
    # plot details
    if PLOT:
        ax.set_xlim(timeE-10, timeCoda + SD + 10)
        outfilename = "{}_{}.png".format(filename, Freq)
        plt.savefig(outfilename)
        plt.close()
    return result


def read_config_file(config_filename):
    config = ConfigParser.SafeConfigParser()
    config.read(config_filename)
    # read main options
    Settings = dict( (k, v) for k, v in config.items("main") )
    Settings['sd'] = int(Settings['sd'])
    # koef
    Settings['koef'] = float(Settings['koef'])
    for k in ("station_lat", "station_lon"): Settings[k] = float(Settings[k])
    # freqs
    Settings["freqs"] = map(float, Settings["freqs"].split())
    # plot (bool)
    Settings["plot"] = True if Settings["plot"].upper() in ("TRUE", "YES") else False
    return Settings


if __name__ == '__main__':
    Settings = read_config_file("coda.conf")
    InputFile = Settings["input_file"]
    if not os.path.exists(InputFile):
        print("Input file to found!")
        sys.exit(1)
    # load input file
    data = np.genfromtxt(InputFile, dtype=None, autostrip=True, names=True,
        loose=True, invalid_raise=False)
    names = data.dtype.names
    # check that all the columns is in the file
    for col in ("TIME_E", "TIME_P", "TIME_S", "FILENAME"):
        if not col in names:
            print("Coulmn {} not found in input file. Check it!".format(col))
            sys.exit(1)
    #===
    WorkBook = xlwt.Workbook()#encoding="utf8")
    Freqs = Settings["freqs"]
    # add sheets with Freq as name
    sheets = [ WorkBook.add_sheet(str(freq)) for freq in Freqs ]
    #=== write headers on every sheet
    if ("station_lat" in Settings.keys() and "station_lon" in Settings.keys()) \
        and ("LAT" in names and "LON" in names):
        has_coords_key = True
        # write column headers: DIST_KM AZIMUTH
        names = list(names) + ["DIST_KM", "DIST_DEGR", "AZIMUTH_AB", "AZIMUTH_BA"]
    else:
        has_coords_key = False
    columns = len(names)
    for sheet in sheets:
        for col, name in enumerate(names):
            sheet.write(0, col, name)
    headers = ["%s_%s" % (w, ch) for ch in ("NS", "EW", "Z") for w in ("P", "S", "C")]
    for sheet in sheets:
        for col, header in enumerate(headers):
            sheet.write(0, col+columns, header)
    #=== main section
    plot = Settings["plot"]
    # view the input file line by line
    for NumLine, line in enumerate(data):
        parts = dict(zip(names, line))
        # analyze parameters
        filename = parts["FILENAME"]
        str_timeE = parts["TIME_E"]
        str_timeP = parts["TIME_P"]
        str_timeS = parts["TIME_S"]
        print " ".join(map(str, line)),
        #== convert time from str to time in seconds (time of Event in seconds since start of day)
        secondsE, secondsP, secondsS = map(calc_seconds_from_time,
            (str_timeE, str_timeP, str_timeS))
        if secondsP < secondsE:
            print("Time P must be after time of Event for {}"),
            continue
        #==
        # for every frequency, calc...
        for NumFreq, Freq in enumerate(Freqs):
            # low and high freq corners
            f1, f2 = 0.5 * Freq, 1.5 * Freq
            #===
            sheet = sheets[NumFreq]
            if has_coords_key:
                # calculate distances, azimuths from the event to the station
                lat, lon = parts["LAT"], parts["LON"]
                # gps2DistAzimuth(lat1, lon1, lat2, lon2)
                dist, azAB, azBA = gps2DistAzimuth(Settings["station_lat"],
                    Settings["station_lon"], lat, lon)
                # locations2degrees(lat1, long1, lat2, long2)
                dist_degr = locations2degrees(Settings["station_lat"],
                    Settings["station_lon"], lat, lon)
                dist /= 1000. # must be in kilometers
                # write distances to sheet
                for col, item in enumerate([dist, dist_degr, azAB, azBA]):
                    sheet.write(NumLine+1, col+columns-4, item)
            # raw data from the directory
            for col, item in enumerate(line):
                sheet.write(NumLine+1, col, item)
            #=== calculations
            try:
                result = main(Freq, f1, f2, filename, secondsE, secondsP,
                    secondsS, Settings)
            except BaseException, e:
                print e,
                result = None
            #=== calculations ended
            if result is not None:
                for col, value in enumerate(result):
                    sheet.write(NumLine+1, col+columns, value)
        print
    # save results
    # filename must consist StationName, algorithm etc
    outfilename = "out_{station}_{algorithm}.xls".format(**Settings)
    try:
        WorkBook.save(outfilename)
    except IOError, e:
        print(e)
