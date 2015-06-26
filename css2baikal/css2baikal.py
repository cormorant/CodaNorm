#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Конвертация данных из формата CSS в Байкал-5.
Всё должно быть в одном файле (нет зависимостей numpy, obspy)
"""
__version__= "0.0.1"
COMPANY_NAME = 'GIN'
APP_NAME = 'css2baikal'

import os
import sys
import struct
import datetime
import argparse
import operator

# edited obspy's CSS
from myobspy.css import isCSS, readCSS
from myobspy.obspycore.stream import Stream

# numpy
import multiarray as ma

#=== Numpy functions
newaxis = None

def asanyarray(a, dtype=None, order=None):
    """ Convert the input to an ndarray, but pass ndarray subclasses through """
    return ma.array(a, dtype, copy=False, order=order, subok=True)


def atleast_3d(*arys):
    """ View inputs as arrays with at least three dimensions """
    res = []
    for ary in arys:
        ary = asanyarray(ary)
        if len(ary.shape) == 0:
            result = ary.reshape(1, 1, 1)
        elif len(ary.shape) == 1:
            result = ary[newaxis,:, newaxis]
        elif len(ary.shape) == 2:
            result = ary[:,:, newaxis]
        else:
            result = ary
        res.append(result)
    if len(res) == 1:
        return res[0]
    else:
        return res


def dstack(tup):
    """ Stack arrays in sequence depth wise (along third axis) """
    return ma.concatenate([atleast_3d(_m) for _m in tup], 2)
#===
#=== Baikal-5 format
# названия каналов
CHANNELS = ("N", "E", "Z")
# полные названия каналов
FULL_CHANNEL_NAMES = ("N-S", "E-W", "Z")
FULL_CHANNEL_NAMES_DICT = dict(zip(CHANNELS, FULL_CHANNEL_NAMES))
#=== Baikal-5 (XX) headers
#= Main header
MainHeaderStruct = '16h16s7d4I'
# имена полей заголовка
MainHeaderNames = ('nkan', 'test', 'vers', 'day', 'month', 'year',
    'satellit', 'valid', 'pri_synhr', 'razr', 'reserv_short1', 'reserv_short2',
    'reserv_short3', 'reserv_short4', 'reserv_short5', 'reserv_short6',
    'station', 'dt', 'to', 'deltas', 'latitude', 'longitude', 'reserv_double1',
    'reserv_double2', 'reserv_long1', 'reserv_long2', 'reserv_long3', 'reserv_long4')
# значения по умолчанию в заголовке
DefaultMainHeaderValues = (0, 0, 50, 1, 1, 2000, 12, 0, 0, 24, 0, 0, 0, 0, 0, 0,
    "st", 0.01, 0., 0., 0., 0., 0., 0., 0, 0, 0, 0)

#= Channel Header
ChannelHeaderStruct = "4h24s24s2d"
# nomer, 3 reserv, name_chan, tip_datch, koef_chan, calcfreq?
ChannelHeaderNames = ('phis_nom', 'reserv1', 'reserv2', 'reserv3',
    'name_chan', 'tip_dat', 'koef_chan', 'reserved')
DefaultChannelHeaderValues = (0, 0, 0, 0, "CH", "t_d", 1.0, 1.0)

FILENAME_FORMAT = "{0.year:04d}-{0.month:02d}-{0.day:02d}_" +\
    "{0.hour:02d}-{0.minute:02d}-{0.second:02d}.xx"
#===


def write_baikal(stream, outdir):
    """ write baikal with metadata """
    # количество каналов/трасс
    nkan = stream.count()
    assert nkan > 0, "Cannot write file with 0 channels!"
    # prepare stream. Here we have "aligned" trace, start/end for traces is same
    starttime = stream[0].stats.starttime
    # where to put final file
    if not os.path.exists(outdir): os.makedirs(outdir)
    outfile = FILENAME_FORMAT.format(starttime)
    outfile = os.path.join(outdir, outfile)
    #=== prepare and write main header
    Dic = dict(zip(MainHeaderNames, DefaultMainHeaderValues))
    Dic["nkan"] = nkan
    # change date (day, month, year)
    date_items = ('day', 'month', 'year')
    getter = operator.attrgetter(*date_items)
    values = getter(starttime)
    for key, value in zip(date_items, values):
        Dic[key] = value
    # calc and change time
    delta = starttime.datetime - datetime.datetime.combine(starttime.date, datetime.time())
    seconds = delta.total_seconds()
    Dic['to'] = seconds
    # station name
    Dic["station"] = str( stream[0].stats.station ) + ".st"
    # save delta
    Dic["dt"] = stream[0].stats.delta
    # write
    values = [Dic[k] for k in MainHeaderNames]
    s_struct = struct.Struct(MainHeaderStruct)
    #==
    # nice output
    s = "Write file from {}.".format(starttime)
    sys.stdout.write("\r" + s)
    sys.stdout.flush()
    #=== baikal-5 file writing starts
    # start writing, do it with Try!
    _f = open(outfile, "wb")
    # write main header
    _f.write(s_struct.pack(*values))
    #= write channel headers
    for num, trace in enumerate(stream):
        # make dictionary with default channel header values
        channelDic = dict(zip(ChannelHeaderNames, DefaultChannelHeaderValues))
        # change number of channel
        channelDic["phis_nom"] = num
        # change name of channel
        channelDic['name_chan'] = FULL_CHANNEL_NAMES_DICT[trace.stats.channel]
        # make koefficient from calibration factor
        channelDic["koef_chan"] = trace.stats.calib
        # write struct channel data
        values = [channelDic[k] for k in ChannelHeaderNames]
        s_struct = struct.Struct(ChannelHeaderStruct)
        _f.write(s_struct.pack(*values))
    # write data
    data = dstack(([trace.data for trace in stream])).flatten()
    _f.write(data.tostring())
    # close
    _f.close()
    return 1


def stream2baikal(stream, args):
    """ запись в формат Байкал-5 из потока с тремя трассами (N, E, Z) по N минут """
    # из минут, сколько это в секундах?
    N = args.minutes * 60
    # начальное время из трех трасс? взять за максимальное из начал
    starttime = max([trace.stats.starttime for trace in stream])
    # конечное время - минимальное из всех
    endtime = min([trace.stats.endtime for trace in stream])
    # дальше цикл, разделяющий на минуты
    end = starttime + N
    # куда записывать результат (в папку со станцией)
    outdir = os.path.join(args.outdir, str(trace.stats.station))
    # before writing, lets change order of channels: N, E, Z
    if (len(stream) == 3) and \
        stream[0].stats.channel.upper() == "E" and \
        stream[1].stats.channel.upper() == "N":
        # change channels
        stream[0], stream[1] = stream[1], stream[0]
    # write in infinite cycle
    while True:
        # save copy, then trim
        newstream = stream.copy()
        newstream.trim(starttime=starttime, endtime=end)#, nearest_sample=False)
        # write it!
        write_baikal(newstream, outdir)
        # time1 now is previous time2; move time2 forward
        starttime = end
        end += N
        # check last
        if end >= endtime:
            # Have to write last samples that left unbounded
            if starttime > endtime: break
            stream.trim(starttime=starttime, endtime=endtime)
            sys.stdout.write("\nWrite from {} to {}".format(starttime, endtime))
            # если последний кучок потока не пустой, записать и его
            if stream:
                write_baikal(stream, outdir)
            break


def process_css(filename, args):
    """ Считывает CSS файл и записывает в ХХ """
    print("Found {}".format(filename))
    # check if it is CSS (header) file
    if not isCSS(filename):
        print("File {} seems not to be a valid CSS file. "
              "Press ENTER to continue anyway or press Ctrl+C to exit.")
        raw_input()
    #===
    # here we have valid CSS file, or continued anyway
    # read CSS header file
    stream = readCSS(filename)
    # fix smth in traces of stream
    for trace in stream:
        #= нужны ли калибровочные коэффициенты. Если нет  - заменим их на 1
        if (not args.savecoef) or (trace.stats.calib == 0.):
            trace.stats.calib = 1.
        # from channel name get 1st 3 elements
        trace.stats.channel = trace.stats.channel[:3]
        #= also rename channels in traces. Should have 1 letter: N, E or Z
        trace.stats.channel = str(trace.stats.channel)
        # check last letter
        if trace.stats.channel[-1].upper() in CHANNELS:
            trace.stats.channel = trace.stats.channel[-1].upper()
        # check 1st letter
        elif trace.stats.channel[0].upper() in CHANNELS:
            trace.stats.channel = trace.stats.channel[0].upper()
        # may be BH1, BH2 (could be also EH1, EH2?)
        elif "1" in trace.stats.channel:
            trace.stats.channel = "N" # 1-й канал это N-S
        elif "2" in trace.stats.channel:
            trace.stats.channel = "E"
        else:
            raise BaseException("What to do with channel name %s? Exiting..." %
                trace.stats.channel)
    # merge, if needed?
    if args.merge:
        # объединить трассы по одной и той же станции и каналу
        try:
            stream.merge(method=1, fill_value='latest')
        except BaseException as e:
            print("Error: {}. Try disabling -s (--savecoef) option.".format(e))
            return
        # sort stream by station, time, channel
        stream.sort(keys=['station', 'channel', 'starttime'])
        # now lets split stream by station
        stations = set([trace.stats.station for trace in stream])
        for station in stations:
            newstream = stream.select(station=station)
            # trim traces on stream, align it
            starttime = max([trace.stats.starttime for trace in newstream])
            endtime = min([trace.stats.endtime for trace in newstream])
            newstream.trim(starttime=starttime, endtime=endtime)
            #===
            # write all traces in stream onto XX files, splitting it, if needed
            stream2baikal(newstream, args)
    else:
        #=== записывать каждую трассу в отдельный файл ХХ, с соотв 1-й точкой...
        # sort stream by time
        stream.sort(keys=['starttime'])
        # куда записывать результат (в папку СТАНЦИЯ/КАНАЛ)
        for trace in stream:
            outdir = os.path.join(args.outdir, str(trace.stats.station),
                str(trace.stats.channel))
            write_baikal(Stream(traces=[trace]), outdir)
        print


def main(args):
    """ main func """
    dirs = args.dirs
    outdir = args.outdir
    if not os.path.exists(outdir): os.makedirs(outdir)
    # can be multiple directories
    for path in args.dirs:
        # проверять, ведь нельзя считывать и записывать одну и ту же папку
        if path == outdir:
            print("Cannot use same path for reading and writing! Skip...")
            continue
        if not os.path.exists(path):
            print("Path %s not found" % path)
            continue
        #===
        # may be it is file
        if os.path.isfile(path):
            process_css(path, args)
        else:
            # find CSS file in path
            for fil in os.listdir(path):
                filename = os.path.join(path, fil)
                # check whether it wfdisc file
                if filename.lower().endswith("wfdisc"):
                    process_css(filename, args)


if __name__ == "__main__":
    #===========================================================================
    # parsers
    parser = argparse.ArgumentParser()
    # общие парсеры
    parser.add_argument('-V', '--version', action='version', 
        version='%(prog)s.' + __version__)
    # конвертировать эти файлы
    parser.add_argument("dirs", nargs='+', help="directories to convert")
    # куда выходные файлы
    parser.add_argument("-o", "--outdir", dest="outdir", default="XX",
        help="path for output data (default is \"XX\")")
    # save calib coef or not
    parser.add_argument("-s", "--savecoef", dest="savecoef", action="store_true",
        help="Save coefficients in CSS file, or set 1.0")
    # merge or not
    parser.add_argument("-m", "--merge", dest="merge", action="store_true",
        help="Merge same channels into one trace with interpolating gaps")
    # по сколько минут писать файлы
    parser.add_argument("-n", "--minutes", dest="minutes", type=int, default=10,
        help="Length of output file (in minutes)")
    #=== endof parsers
    args = parser.parse_args()
    #===========================================================================
    # convert files in arguments, into one file in Baikal-5 format
    try:
        main(args)
    except KeyboardInterrupt:
        print("Interruptted by user!")
        sys.exit(0)
    except (BaseException,) as e:
        print("Error: {}".format(e))
        sys.exit(1)
    else:
        print("\nPress ENTER...")
        raw_input()
