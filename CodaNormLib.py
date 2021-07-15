#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division

"""
Move here code from CodaNorm - for simplicity reason. 
Also moved from Qp module.
"""
import sys
# for Py3 compatibility
if sys.version_info >= (3, ):
    import configparser
else:
    import ConfigParser as configparser


import numpy as np


from statsmodels.robust.robust_linear_model import RLM

from obspy import read_inventory


def linear_fit(y, x, m=None, method='robust', Model=RLM):
    """Code from Qopen: Linear fit between x and y

    :param y,x: data
    :param m: fix slope at specific value
    :param method: one of ('least_squares', 'robust')
    :return: slope a and intercept b of y = ax + b
    """
    if m is None:
        X = np.empty((len(y), 2))
        X[:, 0] = x
        X[:, 1] = 1
        res = Model(y, X).fit()
        return res.params
    else:
        X = np.ones(len(y))
        res = Model(np.array(y) - m * np.array(x), X).fit()
        return m, res.params[0]


def remove_response(stream, station, respfile, pre_filt=[0.001, 0.05, 35, 50], output="VEL"):
    """ get physycal mkm/s units with pre-filter outside of desired frequency range
    pre_filt was [0.001, 0.05, 30, 50] """
    
    inv = read_inventory(respfile)
    stream.remove_response(inventory=inv, output=output, plot=False, pre_filt=pre_filt)
    
    # make physical units mkm/s
    for tr in stream: tr.data *= 1e6
    return stream


def read_settings(config_filename="coda.conf", section="main"):
    config = configparser.ConfigParser()
    config.read(config_filename)
    
    # read main options
    if not config.has_section(section):
        print("No section '{}' in config file {}! Exiting.".format(section,
            config_filename))
        sys.exit(0)
    
    # read all config in dictionary (default type is str)
    Settings = dict( (k, v) for k, v in config.items(section) )
    
    # whats wrong with `channel` option? why always `E`
    Settings['channel'] = config.get(section, "channel")
    
    # read INT value
    for item in ("sd", "sd1", "minsnr", 'corners', "Tcoda"):
        Settings[item] = config.getint(section, item)
    
    # read float values
    for item in ("station_lat", "station_lon", 'koef', "vp", "vs"):
        try:
            Settings[item] = config.getfloat(section, item)
        except configparser.NoOptionError:
            Settings[item] = 0
    # freqs
    Settings["freqs"] = list(map(float, config.get(section, 'freqs').split()))
    
    # read boolean values
    for item in ("plot", 'simulate',  "rms"):#"rotate", "absolute", "envelope",
        Settings[item] = config.getboolean(section, item)
    return Settings


def RMS(array, power=2):
    #return np.sqrt( np.sum(np.power(a.astype(np.float64), power)) / a.size )
    return np.sqrt(np.mean(array ** power))


def calc_signal_noise_ratio(trace, sr, dt, dt2, sd):
    """ calc singal-to-noise ratio 
    Take window for noise at the same time as Time Of Event"""
    
    # slice trace (5 seconds before and 5 after P waves)
    newtrace = trace.slice(dt, dt2)
    
    # take first N seconds of this data for noise, and last N seconds - signal
    # seconds to samples: * sr
    index = int(sd * sr)
    #noise = newtrace.data[100:index+100]# plus one sec
    noise = newtrace.data[int(sr):index+int(sr)]# add one sec
    signal = newtrace.data[-index:]
    
    assert noise.size > 20
    assert noise.size == signal.size
    
    # calc RMS
    return float(RMS(signal) / RMS(noise))


def get_waveforms(stream, network, station, location, channel,
    starttime, endtime, event=None):
    st = stream.select(network=network, station=station, location=location, channel=channel)
    st = st.slice(starttime, endtime)
    return st


def load_stations_from_settingsfile(config_filename="coda.conf", section="stations"):
    """ section 'stations' """
    config = configparser.ConfigParser()
    config.read(config_filename)
    
    # check if exists Stations section
    if not config.has_section(section):
        print("No section '{}' in config file '{}'! Exiting.".format(section,
            config_filename))
        sys.exit(0)
    # read all stations into dictionary (default type is str)
    Stations = dict( (k.upper(), v.split()) for k, v in config.items(section) )
    
    # float coords
    for sta in Stations.keys():
        for N in range(2): Stations[sta][N] = float(Stations[sta][N])
        # must be 4 items for value (STA LAT LON FILENAME DATA)
        errmsg = "Specify 4 items for each STATION: STA LAT LON INPUT_CATALOG DATA_FILENAME"
        assert len(Stations[sta]) == 4, errmsg
    return Stations

#load_stations("stations.dat")
STATIONS = load_stations_from_settingsfile()
