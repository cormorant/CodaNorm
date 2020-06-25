#!/usr/bin/env python
# -*- coding: utf-8 -*-
#from __future__ import division

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
    for item in ("sd", "sd1", "minsnr", 'corners', "r1", "r2"):
        Settings[item] = config.getint(section, item)
    
    # read float values
    for item in ("station_lat", "station_lon", 'koef', "vp", "vs", 'alpha1', 'alpha2'):
        try:
            Settings[item] = config.getfloat(section, item)
        except configparser.NoOptionError:
            Settings[item] = 0
    # freqs
    Settings["freqs"] = list(map(float, config.get(section, 'freqs').split()))
    # plot, rotate? (bool)
    for item in ("plot", 'simulate', "absolute"):#"rotate", 
        Settings[item] = config.getboolean(section, item)
    return Settings


def RMS(array, power=2):
    #return np.sqrt( np.sum(np.power(a.astype(np.float64), power)) / a.size )
    return np.sqrt(np.mean(array ** power))


'''
def get_data_for_window(stream, dt1, dt2):
    """ just get em from data array by seconds. Use slice() method instead! """
    trace = stream[-1]# we use Z channel only
    
    # trim trace from start (P-wave) to end (P-onset + SD)
    newtrace = trace.copy()
    newtrace.trim(dt1, dt2)
    
    # return data
    return newtrace.data
'''

def calc_signal_noise_ratio(stream, dt, dt2, sd):
    """ calc singal-to-noise ratio 
    Take window for noise at the same time as Time Of Event"""
    # use 1 channel only
    trace = stream[-1]
    
    sr = float("%.1f" % trace.stats.sampling_rate)
    
    # slice trace (5 seconds before and 5 after P waves)
    newtrace = trace.slice(dt, dt2)
    
    # take first N seconds of this data for noise, and last N seconds - signal
    # seconds to samples: * sr
    index = int(sd * sr)
    noise = newtrace.data[:index]
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



"""
SD (ZRHB)
SD1=1 for dist=10 -- 23
SD1=2 for dist=15 -- 33
SD1=3 for dist=30 -- 47
SD1=4 for dist=45 -- 58
SD1=5 for dist=53 -- 75
SD1=6 for dist=69 -- 89
SD1=7 for dist=100-- 98
SD1=8 for dist=102-- 99
SD1=9 for dist=119-- 130
SD1=10 for dist=136- 149
SD1=11 for dist=153- 165
SD1=12 for dist=156
SD1=13 for dist=192
SD1=14 for dist=188- 206
SD1=15 for dist=197- 209
SD1=16 for dist=217- 237
"""
