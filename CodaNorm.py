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
Программа расчета коды сейсмограммы по всем (?) каналам.
Входной формат данных - GSE2.

Расчет коды осуществляется по 2-м алгоритмам:
    - алгоритм нормализации.
    - по амплитудным спектрам коды, волн P, S.
"""
APP_NAME = "CodaNorm"
__version__="0.1"
COMPANY_NAME = 'GIN SB RAS'

import os
import sys
from obspy.gse2.core import isGSE2, readGSE2
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
except ImportError, e:
    print e
    matplotlib_is_imported = False
else:
    matplotlib_is_imported = True
    COLORS = cycle(("g", "c", "r"))


def calc_spectrum(y, dt):
    """ calc spectrum """
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
    """ вычисляет время в секундах с начала дня из текстового времени """
    hour, minute, second = map(float, _timeStr.split(":"))
    return hour * 3600 + minute * 60 + second


def calc_power(P, S, C):
    """ 
    На входе - спектры в окнах отфильтрованных значений (P, S, кода).
    Среднеквадратическая оценка спектра
    """
    return [np.sqrt( np.sum(np.power(a, 2)) / a.size ) for a in (P, S, C)]


def main(Freq, f1, f2, filename, secondsE, secondsP, secondsS, Settings, rotate=None):
    """ расчет параметров коды """
    # plot or not
    PLOT = Settings["plot"]
    if not matplotlib_is_imported: PLOT = False
    # LENGTH OF WINDOWS and KOEF
    SD = Settings["sd"]
    KOEF = Settings["koef"]
    # открывать исходный файл
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
    # будем считать спектры по окну P, S, коде
    for ax, trace in zip(axes, stream.traces):
        sr = trace.stats.sampling_rate # sampling rate, 100 Hz normally
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
        # за 0 берем время начала файла
        secondsStart = calc_seconds_from_time(trace.stats.starttime.time.strftime("%H:%M:%S.%f"))
        # assert check, если файл начинается после события - что делать?
        if secondsStart > secondsE:
            print("File {} starts after Event!"),
            return
        # уменьшим времена в секундах относительно начала файла
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
        #=== все вычисления по данной компоненте
        #= window for P
        Pindex1 = int( np.ceil(timeP * sr) )
        Pindex2 = int( Pindex1 + SD * sr )
        Pwindow = y[Pindex1:Pindex2] # P
        #= window for S
        Sindex1 = int( np.ceil(timeS*sr) )
        Sindex2 = int( Sindex1 + SD * sr )
        Swindow = y[Sindex1:Sindex2] # S
        # C - кода
        #TODO: add an option to change this multiplier (originally 2.5)
        coda_index1 = int(np.ceil(KOEF * (timeS - timeE)) * sr + np.ceil(timeE * sr))
        coda_index2 = int(coda_index1 + SD * sr)
        # вырезаный по индексам массив с кодой
        coda_window = y[coda_index1:coda_index2]
        if not coda_window.size:
            print("Not enough time for coda to cut!"),
            return
        #=== spectrums or just normalize
        # считать средне-квадр. значения по спектру
        if "SPECTR" in Settings["algorithm"].upper():
            freqs, valuesP = calc_spectrum(Pwindow, trace.stats.delta) # P
            freqs, valuesS = calc_spectrum(Swindow, trace.stats.delta) # S
            freqs, valuesC = calc_spectrum(coda_window, trace.stats.delta) # C
        # или по исходным данным
        else:
            valuesP, valuesS, valuesC = Pwindow, Swindow, coda_window
        # calc result on windows
        result += calc_power(valuesP, valuesS, valuesC)
    # plot details
    if PLOT:
        ax.set_xlim(timeE-10, timeCoda + SD + 10)
        #plt.show()
        outfilename = "{}_{}.png".format(filename, Freq)
        plt.savefig(outfilename)
        plt.close()
    return result


def read_config_file(config_filename):
    config = ConfigParser.SafeConfigParser()
    config.read(config_filename)
    # read main options
    Settings = dict( (k, v) for k, v in config.items("main") )
    # длина окна - целое число
    Settings['sd'] = int(Settings['sd'])
    # koef
    Settings['koef'] = float(Settings['koef'])
    # координаты - веществ число
    for k in ("station_lat", "station_lon"): Settings[k] = float(Settings[k])
    # freqs
    Settings["freqs"] = map(float, Settings["freqs"].split())
    # plot (bool)
    Settings["plot"] = True if Settings["plot"].upper() in ("TRUE", "YES") else False
    return Settings


if __name__ == '__main__':
    #===
    # Для каждой записи получить:
    # времена события (задается во входном каталоге)
    # - время вступления волны P для данной станции
    # - время вступления волны S для данной станции
    # - имя файла для обработки, откуда считывать данные
    #=== считывание файла настроек
    #try:
    Settings = read_config_file("coda.conf")
    #except NoSectionError
    InputFile = Settings["input_file"]
    if not os.path.exists(InputFile):
        print("Input file to found!")
        sys.exit(1)
    # загрузить входной файл (guess types in file, hoping for the best)
    data = np.genfromtxt(InputFile, dtype=None, autostrip=True, names=True,
        loose=True, invalid_raise=False)
    names = data.dtype.names
    # проверить что все нужные столбцы в файле есть
    for col in ("TIME_E", "TIME_P", "TIME_S", "FILENAME"):
        if not col in names:
            print("Coulmn {} not found in input file. Check it!".format(col))
            sys.exit(1)
    #===
    # собирать значения и сохранять в файл Эксель
    WorkBook = xlwt.Workbook()#encoding="utf8")
    Freqs = Settings["freqs"]
    # add sheets with Freq as name
    sheets = [ WorkBook.add_sheet(str(freq)) for freq in Freqs ]
    #=== write headers on every sheet
    # если есть настройки с координатами, и во входном файле есть коорд
    if ("station_lat" in Settings.keys() and "station_lon" in Settings.keys()) \
        and ("LAT" in names and "LON" in names):
        has_coords_key = True
        # write column headers: DIST_KM AZIMUTH
        names = list(names) + ["DIST_KM", "DIST_DEGR", "AZIMUTH_AB", "AZIMUTH_BA"]
    else:
        has_coords_key = False
    # сколько столбцов занято под исходные данные = кол-во names
    columns = len(names)
    # заголовки на каждом листе
    for sheet in sheets:
        for col, name in enumerate(names):
            sheet.write(0, col, name)
    # записать остальные заголовки для каналов и названий окон P, S, C
    # три канала, три окна: P_N, S_N, C_N
    headers = ["%s_%s" % (w, ch) for ch in ("NS", "EW", "Z") for w in ("P", "S", "C")]
    for sheet in sheets:
        for col, header in enumerate(headers):
            sheet.write(0, col+columns, header)
    #=== main section
    plot = Settings["plot"]
    # просматривать входной файл построчно
    for NumLine, line in enumerate(data):
        # обрабатывать данные по событию
        parts = dict(zip(names, line))
        # разобрать параметры
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
        # для всех частот считать...
        for NumFreq, Freq in enumerate(Freqs):
            # подсчитать границы от центральной частоты, low and high freq corners
            f1, f2 = 0.5 * Freq, 1.5 * Freq
            #===
            # подготовить, куда писать результаты
            sheet = sheets[NumFreq]
            if has_coords_key:
                # считать расстояние, азимуты от события до станции (или наоборот?)
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
            # исходные данные из каталога
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
            # выводить значения по каналам, результат
            if result is not None:
                #print " ".join(map(str, result)),
                for col, value in enumerate(result):
                    sheet.write(NumLine+1, col+columns, value)
        print
    # save results
    # filename must consist StationName, alprithm
    outfilename = "out_{station}_{algorithm}.xls".format(**Settings)
    try:
        WorkBook.save(outfilename)
        #print
    except IOError, e:
        print(e)
