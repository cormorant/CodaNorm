#coding: utf-8
"""
Convert baikal file to seisan format using obspy
"""
import os
import sys
import datetime
from obspy.core import Trace, Stream, UTCDateTime
import numpy as np

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


if __name__ == "__main__":
    # в указанной папке конвертировать все файлы
    # сохраняя результат в соотв папке gse2
    Dir = "data"
    for root, dirnames, filenames in os.walk(Dir):
        for filename in filenames:
            filename = os.path.join(root, filename)
            print("Converting file %s" % filename)
            try: write_seisan(filename)
            except: print("Error writing file %s" % filename)
