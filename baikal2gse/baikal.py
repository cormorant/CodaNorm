#!/usr/bin/env python
# coding: utf-8
from __future__ import division
"""
Description of XX (Baikal-5) format
"""
APP_NAME = "BaikalLib"
__version__="0.0.1.2"
COMPANY_NAME = 'GIN'

import os
import sys
import struct
import numpy as np

# short names
MainHeaderNames = ('nkan', 'test', 'vers', 'day', 'month', 'year',
    'satellit', 'valid', 'pri_synhr', 'razr', 'reserv_short1', 'reserv_short2',
    'reserv_short3', 'reserv_short4', 'reserv_short5', 'reserv_short6',
    'station', 'dt', 'to', 'deltas', 'latitude', 'longitude')
MainHeaderTypeStruct = '16h16s5d'

# channel names
ChannelHeaderNames = ('phis_nom', 'reserv1', 'reserv2', 'reserv3',
    'name_chan', 'tip_dat', 'koef_chan', 'reserved')
ChannelHeaderTypeStruct = '4h24s24s2d'

# start of channel headers
CHANNEL_HEADER_START_OFFSET = 120


class BaikalFile():
    """ XX (Baikal-5) format """
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename, "rb") as _f:
            if self._is_baikal(_f):
                self.valid = True
                # read main header
                self.MainHeader = self._readMainHeader(_f)
                # read channel header
                self.ChannelHeaders = self._readChannelHeaders(_f)
                # read data array
                self.traces = self._readData(_f)
            else:
                self.valid = False

    def _is_baikal(self, _f):
        """ is valid baikal file or not """
        _f.seek(0)
        try:
            nkan = struct.unpack("h", _f.read(2))[0]
        except struct.error:
            print("struct error at %s" % self.filename)
            return
        if not nkan in range(1,7): return
        return True

    def stripnulls(self, s):
        """ clean string buffer """
        s = s.strip()
        for sym in ("\00", "\01", ".st"): s = s.strip(sym)
        return s

    def _readMainHeader(self, _f):
        """ read main header """
        _f.seek(0)
        size = struct.calcsize(MainHeaderTypeStruct)
        data = struct.unpack(MainHeaderTypeStruct, _f.read(size))
        header = dict(zip(MainHeaderNames, data))
        # fix station name
        header["station"] = self.stripnulls(header["station"])
        # fix invalid year, if possible
        if header["year"] < 1900: header["year"] += 2000
        return header

    def _readChannelHeaders(self, _f):
        """ read channel header """
        start = CHANNEL_HEADER_START_OFFSET
        nkan = self.MainHeader['nkan']
        size = struct.calcsize(ChannelHeaderTypeStruct)
        headers = []
        _f.seek(start)
        for kan in range(nkan):
            data = struct.unpack(ChannelHeaderTypeStruct, _f.read(size))
            result = dict(zip(ChannelHeaderNames, data))
            result["name_chan"] = self.stripnulls( result["name_chan"] )
            result['tip_dat'] = self.stripnulls( result['tip_dat'] )
            headers += [ result ]
        return headers
    
    def _readData(self, _f):
        """ read data """
        nkan = self.MainHeader['nkan']
        razr = self.MainHeader['razr']
        razm = 2 if razr==16 else 4
        # type
        typ = "h" if razr==16 else "i"
        # dtype
        dtyp = np.int16 if razr==16 else np.int32
        # start offset
        offset = 120 + nkan * 72
        # load and read
        _f.seek(offset)
        data = np.fromstring(_f.read(), dtype=dtyp)
        # strip array
        while data.shape[0] % nkan != 0:
            data = data[:-1]
        # demultiplexing
        data = data.reshape(int(data.shape[0]/nkan), nkan).T
        return data

    def get_time(self, to):
        """ Return time (HH, MM, SS) from seconds since the beginning of the day """
        hours, remainder = divmod(to, 3600)
        minutes, seconds = divmod(remainder, 60)
        return hours, minutes, seconds
