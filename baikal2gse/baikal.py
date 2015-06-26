#!/usr/bin/env python
# coding: utf-8
from __future__ import division
"""
Описание формата Байкал-5
"""
APP_NAME = "BaikalLib"
__version__="0.0.1.2"
COMPANY_NAME = 'GIN'

import os
import sys
import struct
import numpy as np

# Краткое описание типов данных:
# имена
MainHeaderNames = ('nkan', 'test', 'vers', 'day', 'month', 'year',
    'satellit', 'valid', 'pri_synhr', 'razr', 'reserv_short1', 'reserv_short2',
    'reserv_short3', 'reserv_short4', 'reserv_short5', 'reserv_short6',
    'station', 'dt', 'to', 'deltas', 'latitude', 'longitude')
# заголовок файла
MainHeaderTypeStruct = '16h16s5d'

# имена
ChannelHeaderNames = ('phis_nom', 'reserv1', 'reserv2', 'reserv3',
    'name_chan', 'tip_dat', 'koef_chan', 'reserved')
# заголовок канала:
ChannelHeaderTypeStruct = '4h24s24s2d'

# адрес начала структур заголовка каналов (и размер главного заголовка)
CHANNEL_HEADER_START_OFFSET = 120


class BaikalFile():
    """ Описание формата Байкал-5 """
    def __init__(self, filename):
        self.filename = filename
        with open(self.filename, "rb") as _f:
            # если файл является файлом формата Байкал - работаем с ним
            if self._is_baikal(_f):
                self.valid = True
                # считаем заголовок
                self.MainHeader = self._readMainHeader(_f)
                # заголовки каналов
                self.ChannelHeaders = self._readChannelHeaders(_f)
                # читать массивы с данными
                self.traces = self._readData(_f)
            else:
                self.valid = False

    def _is_baikal(self, _f):
        """ является ли файлом формата Байкал """
        _f.seek(0)
        # количество каналов
        try:
            nkan = struct.unpack("h", _f.read(2))[0]
        except struct.error:
            print("struct error at %s" % self.filename)
            return
        # должно быть вразумительное число каналов
        if not nkan in range(1,7): return
        #TODO: проверка правильности даты, разрядности и тд
        # если сюда дошло - все проверки выполнены
        return True

    def stripnulls(self, s):
        """ очищает строку от символов пропуска и нулевых символов """
        #TODO: удалять цифры из строки станции (re.compile...)
        s = s.strip()
        for sym in ("\00", "\01", ".st"): s = s.strip(sym)
        return s

    def _readMainHeader(self, _f):
        """ считывание заголовка файла """
        # читать все записи из структуры или только обязательные
        _f.seek(0)
        size = struct.calcsize(MainHeaderTypeStruct)#120 if all
        data = struct.unpack(MainHeaderTypeStruct, _f.read(size))
        header = dict(zip(MainHeaderNames, data))
        # поправим станцию
        header["station"] = self.stripnulls(header["station"])
        # неправильный год кое-где
        if header["year"] < 1900: header["year"] += 2000
        return header

    def _readChannelHeaders(self, _f):
        """ считывание заголовков каналов, структура CHANNEL_HEADER """
        start = CHANNEL_HEADER_START_OFFSET
        nkan = self.MainHeader['nkan']
        size = struct.calcsize(ChannelHeaderTypeStruct)
        headers = []
        _f.seek(start)
        for kan in range(nkan):
            # считывание очередного канала
            data = struct.unpack(ChannelHeaderTypeStruct, _f.read(size))
            result = dict(zip(ChannelHeaderNames, data))
            result["name_chan"] = self.stripnulls( result["name_chan"] )
            result['tip_dat'] = self.stripnulls( result['tip_dat'] )
            headers += [ result ]
        return headers
    
    def _readData(self, _f):
        """ считывание данных """
        nkan = self.MainHeader['nkan']
        razr = self.MainHeader['razr']
        # размер одного замера
        razm = 2 if razr==16 else 4
        # тип
        typ = "h" if razr==16 else "i"
        # dtype
        dtyp = np.int16 if razr==16 else np.int32
        # где начинать считывать данные (336)
        offset = 120 + nkan * 72
        # load&read
        _f.seek(offset)
        data = np.fromstring(_f.read(), dtype=dtyp)
        # обрезать массив с конца пока он не делится на 3
        while data.shape[0] % nkan != 0:
            data = data[:-1]
        # демультиплексируем
        data = data.reshape(int(data.shape[0]/nkan), nkan).T
        return data

    def get_time(self, to):
        """ Возвращаем вычисленное время (Ч, М, С) из числа секунд с начала суток """
        hours, remainder = divmod(to, 3600)
        minutes, seconds = divmod(remainder, 60)
        # вернём время
        return hours, minutes, seconds
