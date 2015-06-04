#coding: utf-8
"""
Convert baikal file to seisan format
"""
__version__="0.0.1"
COMPANY_NAME = 'GIN'
APP_NAME = 'baikal2gse'

import os
import sys
import datetime
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np
import argparse

from baikal import BaikalFile

import re
_split = re.compile(r'[\0%s]' %
    re.escape(''.join([os.path.sep, os.path.altsep or ''])))

def secure_filename(path):
    return _split.sub('', path)


DEFAULT_STATS = {
    'network': "",
    'location': "LOC",
    "calib": 1.0,
}

CHANNELS = ("N", "E", "Z")


def write_seisan(filename):
    """ writes seisan file from baikal one """
    bf = BaikalFile(filename)
    if not bf.valid:
        print("Invalid file {}".format(filename))
        return
    header = bf.MainHeader
    # datetime
    date = datetime.datetime(header["year"], header["month"], header["day"])
    delta = datetime.timedelta(seconds=header["to"])
    dt = date + delta
    _time = dt.time() # time
    # make utc datetime
    utcdatetime = UTCDateTime(date.year, date.month, date.day,
        _time.hour, _time.minute, _time.second, _time.microsecond, precision=3)
    # названия каналов в выходном файле (N, E, Z)
    bf.traces = bf.traces.astype(np.int32)
    #! Берём первые три (3) канала. Грубые каналы отбрасываются
    bf.traces = bf.traces[:3]
    traces = []
    for channel, data in zip(CHANNELS, bf.traces):
        # подготовить заголовок
        stats = DEFAULT_STATS.copy()
        stats.update({
            "station": header['station'].upper()[:3],
            'channel': channel,
            'sampling_rate': int( 1./header["dt"] ),
            "delta": header["dt"],
            "npts": data.size,#shape[0]
            'starttime': utcdatetime,
        })
        # создать трассу
        trace = Trace(data=data, header=stats)
        # объединять все в одну трассу
        traces.append(trace)
    # create Stream
    stream = Stream(traces)
    # write seisan
    # path
    path = os.path.join("gse2", stats["station"], str(stats['starttime'].year))
    path = secure_filename(path)
    if not os.path.exists(path): os.makedirs(path)
    #=== filename (using format)
    # date
    name = "{year:04}-{month:02}-{day:02}".format(**header)
    # time
    name += "-{t.hour:02}-{t.minute:02}".format(t=stats['starttime'])
    # + station name + Day_of_Year
    name += "{0}__{1:03}".format(stats["station"], stats['starttime'].timetuple().tm_yday)
    print 'writing', name
    stream.write(os.path.join(path, name), format='GSE2')


def main(args):
    """ convert files (or directories) to gse2 format """
    # parse arguments
    dirs = args.dirs
    outdir = args.outdir
    if not os.path.exists(outdir): os.makedirs(outdir)
    '''
    # times
    starttime = args.start
    endtime = args.end
    # convert to UTCDateTime
    starttime = str2UTCDateTime(starttime)
    # decrease starttime in 1 minute?
    endtime = str2UTCDateTime(endtime)
    '''
    # create Stream
    stream = Stream()
    # read mseed files
    for filename in files:
        stream += readMSEED(filename, format="MSEED", starttime=starttime, endtime=endtime)
    # may be there no times: then align stream!
    starttime = max([trace.stats.starttime for trace in stream])
    endtime = min([trace.stats.endtime for trace in stream])
    stream.trim(starttime=starttime, endtime=endtime)
    # write it!
    print "Writing {network}_{station}_{location}({starttime}-{endtime})".format(**stream[0].stats)
    result = write_baikal(stream, args)


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
    parser.add_argument("-o", "--outdir", dest="outdir", default="gse2",
        help="path for output data (default is \"gse2\")")
    # start from
    #parser.add_argument("-s", "--start", dest="start",
    #    help="start time (format \"YYYY-mm-ddTHH:MM:SS\")"
    #    " like 2015-05-28T07:12:56")
    # end time
    #parser.add_argument("-e", "--end", dest="end",
    #    help="end time (format YYYY-mm-ddTHH:MM:SS)")
    #== end of parsing parameters
    args = parser.parse_args()
    print args
    sys.exit(0)
    #===========================================================================
    #
    #TODO: проверять, ведь нельзя считывать и записывать одну и ту же папку
    #
    #=== start job
    # convert files in arguments, into one file in Baikal-5 format
    try:
        main(args)
    except (ValueError,) as e:
        print("Error: {}".format(e))
        sys.exit(1)

