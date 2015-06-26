# -*- coding: utf-8 -*-
"""
Module for handling ObsPy Trace objects.
"""
from __future__ import division
import functools
import math
import warnings
from copy import copy, deepcopy

import numpy as np

#import compatibility
from mygse.compatibility import round_away, getcallargs

from mygse.obspycore.utcdatetime import UTCDateTime
from mygse.obspycore.attribdict import AttribDict


class Stats(AttribDict):
    """ A container for additional header information of a ObsPy Trace object """
    readonly = ['endtime']
    defaults = {
        'sampling_rate': 1.0,
        'delta': 1.0,
        'starttime': UTCDateTime(0),
        'endtime': UTCDateTime(0),
        'npts': 0,
        'calib': 1.0,
        'network': '',
        'station': '',
        'location': '',
        'channel': '',
    }

    def __init__(self, header={}):
        """
        """
        super(Stats, self).__init__(header)

    def __setitem__(self, key, value):
        """
        """
        # keys which need to refresh derived values
        if key in ['delta', 'sampling_rate', 'starttime', 'npts']:
            # ensure correct data type
            if key == 'delta':
                key = 'sampling_rate'
                value = 1.0 / float(value)
            elif key == 'sampling_rate':
                value = float(value)
            elif key == 'starttime':
                value = UTCDateTime(value)
            elif key == 'npts':
                value = int(value)
            # set current key
            super(Stats, self).__setitem__(key, value)
            # set derived value: delta
            try:
                delta = 1.0 / float(self.sampling_rate)
            except ZeroDivisionError:
                delta = 0
            self.__dict__['delta'] = delta
            # set derived value: endtime
            if self.npts == 0:
                timediff = 0
            else:
                timediff = (self.npts - 1) * delta
            self.__dict__['endtime'] = self.starttime + timediff
            return
        # prevent a calibration factor of 0
        if key == 'calib' and value == 0:
            msg = 'Calibration factor set to 0.0!'
            warnings.warn(msg, UserWarning)
        # all other keys
        if isinstance(value, dict):
            super(Stats, self).__setitem__(key, AttribDict(value))
        else:
            super(Stats, self).__setitem__(key, value)

    __setattr__ = __setitem__

    def __str__(self):
        """
        Return better readable string representation of Stats object.
        """
        priorized_keys = ['network', 'station', 'location', 'channel',
                          'starttime', 'endtime', 'sampling_rate', 'delta',
                          'npts', 'calib']
        return self._pretty_str(priorized_keys)

    def _repr_pretty_(self, p, cycle):
        p.text(str(self))


def _add_processing_info(func):
    """
    This is a decorator that attaches information about a processing call as a
    string to the Trace.stats.processing list.
    """
    @functools.wraps(func)
    def new_func(*args, **kwargs):
        callargs = getcallargs(func, *args, **kwargs)
        callargs.pop("self")
        kwargs_ = callargs.pop("kwargs", {})
        info = "ObsPy {0.0}: function"
        arguments = []
        arguments += \
            ["%s=%s" % (k, v) if not isinstance(v, str) else
             "%s='%s'" % (k, v) for k, v in callargs.items()]
        arguments += \
            ["%s=%s" % (k, v) if not isinstance(v, str) else
             "%s='%s'" % (k, v) for k, v in kwargs_.items()]
        arguments.sort()
        info = info# % "::".join(arguments)
        self = args[0]
        result = func(*args, **kwargs)
        # Attach after executing the function to avoid having it attached
        # while the operation failed.
        self._addProcessingInfo(info)
        return result

    new_func.__name__ = func.__name__
    new_func.__doc__ = func.__doc__
    new_func.__dict__.update(func.__dict__)
    return new_func


class Trace(object):
    """ An object containing data of a continuous series, such as a seismic trace """

    def __init__(self, data=np.array([]), header=None):
        # make sure Trace gets initialized with suitable ndarray as self.data
        # otherwise we could end up with e.g. a list object in self.data
        _data_sanity_checks(data)
        # set some defaults if not set yet
        if header is None:
            # Default values: For detail see
            # http://www.obspy.org/wiki/\
            # KnownIssues#DefaultParameterValuesinPython
            header = {}
        header.setdefault('npts', len(data))
        self.stats = Stats(header)
        # set data without changing npts in stats object (for headonly option)
        super(Trace, self).__setattr__('data', data)

    @property
    def meta(self):
        return self.stats

    @meta.setter
    def meta(self, value):
        self.stats = value

    def __eq__(self, other):
        """
        Implements rich comparison of Trace objects for "==" operator.

        Traces are the same, if both their data and stats are the same.
        """
        # check if other object is a Trace
        if not isinstance(other, Trace):
            return False
        # comparison of Stats objects is supported by underlying AttribDict
        if not self.stats == other.stats:
            return False
        # comparison of ndarrays is supported by NumPy
        if not np.array_equal(self, other):
            return False

        return True

    def __ne__(self, other):
        """
        Implements rich comparison of Trace objects for "!=" operator.

        Calls __eq__() and returns the opposite.
        """
        return not self.__eq__(other)

    def __len__(self):
        return len(self.data)

    count = __len__

    def __setattr__(self, key, value):
        """
        __setattr__ method of Trace object.
        """
        # any change in Trace.data will dynamically set Trace.stats.npts
        if key == 'data':
            _data_sanity_checks(value)
            self.stats.npts = len(value)
        return super(Trace, self).__setattr__(key, value)

    def __getitem__(self, index):
        """
        __getitem__ method of Trace object.

        :rtype: list
        :return: List of data points
        """
        return self.data[index]

    def getId(self):
        """ Return a SEED compatible identifier of the trace """
        out = "%(network)s.%(station)s.%(location)s.%(channel)s"
        return out % (self.stats)

    id = property(getId)

    def _ltrim(self, starttime, pad=False, nearest_sample=True,
               fill_value=None):
        """
        Cut current trace to given start time. For more info see
        :meth:`~obspy.core.trace.Trace.trim`.

        .. rubric:: Example

        >>> tr = Trace(data=np.arange(0, 10))
        >>> tr.stats.delta = 1.0
        >>> tr._ltrim(tr.stats.starttime + 8)  # doctest: +ELLIPSIS
        <...Trace object at 0x...>
        >>> tr.data
        array([8, 9])
        >>> tr.stats.starttime
        UTCDateTime(1970, 1, 1, 0, 0, 8)
        """
        org_dtype = self.data.dtype
        if isinstance(starttime, float) or isinstance(starttime, int):
            starttime = UTCDateTime(self.stats.starttime) + starttime
        elif not isinstance(starttime, UTCDateTime):
            raise TypeError
        # check if in boundary
        if nearest_sample:
            delta = round_away((starttime - self.stats.starttime) * self.stats.sampling_rate)
            # due to rounding and npts starttime must always be right of
            # self.stats.starttime, rtrim relies on it
            if delta < 0 and pad:
                npts = abs(delta) + 10  # use this as a start
                newstarttime = self.stats.starttime - npts / \
                    float(self.stats.sampling_rate)
                newdelta = round_away((starttime - newstarttime) * self.stats.sampling_rate)
                delta = newdelta - npts
            delta = int(delta)
        else:
            delta = int(math.floor(round((self.stats.starttime - starttime) *
                                   self.stats.sampling_rate, 7))) * -1
        # Adjust starttime only if delta is greater than zero or if the values
        # are padded with masked arrays.
        if delta > 0 or pad:
            self.stats.starttime += delta * self.stats.delta
        if delta == 0 or (delta < 0 and not pad):
            return self
        elif delta < 0 and pad:
            try:
                gap = createEmptyDataChunk(abs(delta), self.data.dtype,
                                           fill_value)
            except ValueError:
                # createEmptyDataChunk returns negative ValueError ?? for
                # too large number of points, e.g. 189336539799
                raise Exception("Time offset between starttime and "
                                "trace.starttime too large")
            self.data = np.ma.concatenate((gap, self.data))
            return self
        elif starttime > self.stats.endtime:
            self.data = np.empty(0, dtype=org_dtype)
            return self
        elif delta > 0:
            try:
                self.data = self.data[delta:]
            except IndexError:
                # a huge numbers for delta raises an IndexError
                # here we just create empty array with same dtype
                self.data = np.empty(0, dtype=org_dtype)
        return self

    def _rtrim(self, endtime, pad=False, nearest_sample=True, fill_value=None):
        """
        Cut current trace to given end time. For more info see
        :meth:`~obspy.core.trace.Trace.trim`.

        .. rubric:: Example

        >>> tr = Trace(data=np.arange(0, 10))
        >>> tr.stats.delta = 1.0
        >>> tr._rtrim(tr.stats.starttime + 2)  # doctest: +ELLIPSIS
        <...Trace object at 0x...>
        >>> tr.data
        array([0, 1, 2])
        >>> tr.stats.endtime
        UTCDateTime(1970, 1, 1, 0, 0, 2)
        """
        org_dtype = self.data.dtype
        if isinstance(endtime, float) or isinstance(endtime, int):
            endtime = UTCDateTime(self.stats.endtime) - endtime
        elif not isinstance(endtime, UTCDateTime):
            raise TypeError
        # check if in boundary
        if nearest_sample:
            delta = round_away((endtime - self.stats.starttime) *
                self.stats.sampling_rate) - self.stats.npts + 1
            delta = int(delta)
        else:
            # solution for #127, however some tests need to be changed
            # delta = -1*int(math.floor(compatibility.round_away(
            #     (self.stats.endtime - endtime) * \
            #     self.stats.sampling_rate, 7)))
            delta = int(math.floor(round((endtime - self.stats.endtime) *
                                   self.stats.sampling_rate, 7)))
        if delta == 0 or (delta > 0 and not pad):
            return self
        if delta > 0 and pad:
            try:
                gap = createEmptyDataChunk(delta, self.data.dtype, fill_value)
            except ValueError:
                # createEmptyDataChunk returns negative ValueError ?? for
                # too large number of points, e.g. 189336539799
                raise Exception("Time offset between starttime and " +
                                "trace.starttime too large")
            self.data = np.ma.concatenate((self.data, gap))
            return self
        elif endtime < self.stats.starttime:
            self.stats.starttime = self.stats.endtime + \
                delta * self.stats.delta
            self.data = np.empty(0, dtype=org_dtype)
            return self
        # cut from right
        delta = abs(delta)
        total = len(self.data) - delta
        if endtime == self.stats.starttime:
            total = 1
        self.data = self.data[:total]
        return self

    @_add_processing_info
    def trim(self, starttime=None, endtime=None, pad=False,
             nearest_sample=True, fill_value=None):
        """ Cut current trace to given start and end time """
        # check time order and swap eventually
        if starttime and endtime and starttime > endtime:
            raise ValueError("startime is larger than endtime")
        # cut it
        if starttime:
            self._ltrim(starttime, pad, nearest_sample=nearest_sample,
                        fill_value=fill_value)
        if endtime:
            self._rtrim(endtime, pad, nearest_sample=nearest_sample,
                        fill_value=fill_value)
        # if pad=True and fill_value is given convert to NumPy ndarray
        if pad is True and fill_value is not None:
            try:
                self.data = self.data.filled()
            except AttributeError:
                # numpy.ndarray object has no attribute 'filled' - ignoring
                pass
        return self

    def max(self):
        """ Returns the value of the absolute maximum amplitude in the trace """
        value = self.data.max()
        _min = self.data.min()
        if abs(_min) > abs(value):
            value = _min
        return value

    def copy(self):
        """ Returns a deepcopy of the trace """
        return deepcopy(self)

    def _addProcessingInfo(self, info):
        """
        Add the given informational string to the `processing` field in the
        trace's :class:`~obspy.core.trace.Stats` object.
        """
        proc = self.stats.setdefault('processing', [])
        proc.append(info)

    def times(self):
        """ For convenient plotting compute a NumPy array of seconds since
        starttime corresponding to the samples in Trace """
        timeArray = np.arange(self.stats.npts)
        timeArray = timeArray / self.stats.sampling_rate
        # Check if the data is a ma.maskedarray
        if isinstance(self.data, np.ma.masked_array):
            timeArray = np.ma.array(timeArray, mask=self.data.mask)
        return timeArray


def _data_sanity_checks(value):
    """
    Check if a given input is suitable to be used for Trace.data. Raises the
    corresponding exception if it is not, otherwise silently passes.
    """
    if not isinstance(value, np.ndarray):
        msg = "Trace.data must be a NumPy array."
        raise ValueError(msg)
    if value.ndim != 1:
        msg = ("NumPy array for Trace.data has bad shape ('%s'). Only 1-d "
               "arrays are allowed for initialization.") % str(value.shape)
        raise ValueError(msg)
