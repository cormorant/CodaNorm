#coding: utf-8
"""
Ищем события в бд-acess.
"""
APP_NAME = "find_eq"
__version__="0.0.1.2"
COMPANY_NAME = 'GIN'

import os
import sys
import pyodbc
import datetime
import numpy as np

from BaikalFile import BaikalFile


# откуда считать секунды в  UTevent
DT0 = datetime.datetime(2000, 1, 1)

EQ_FILE = "catalog.dat"
# сколько брать +-секунд? расстояние до 300 км, скорость 4км/с = ~50 c
INTERVAL = 110

# database #1 (2001-2005 years)
DB_FILE1 = "D:/Work/seis/seisobr/PrnBase01_05.mdb"
# database #2 (2006-2011 years)
DB_FILE2 = "D:/Work/seis/seisobr/PrnBase2006-2011.mdb"

#=== sql queries
SELECT_PrnsDir = '''\
SELECT idDir, Path, UTevent, DateE, TimeE, Energy
FROM PrnsDir
WHERE DateE=? AND UTevent BETWEEN ? AND ?\
'''

SELECT_Prns = '''\
SELECT idPrn, idDir, FName, seisFile, K, seismgrStName
FROM Prns
WHERE idDir=?\
'''

SELECT_WAVES = """\
SELECT NameWave, TimeWave, komp, Z
FROM PrnsWaves
WHERE idPrn=?
AND NameWave LIKE '__m'\
"""

STATIONS = (
    "HRM", 
    #"UUD", 
    #"TRT",
    #"ZRH",
    #"MXM",
    #"FFN",
    #"KEL",
    #"VBR",
)


def setup_db_conn(mdbfile):
    # conAcc = pyodbc.connect(
    #          'DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=Northwind.accdb')
    #conn_str = "DRIVER={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s; Provider=MSDASQL;" % mdbfile
    conn_str = "Driver={Microsoft Access Driver (*.mdb, *.accdb)};DBQ=%s;" % mdbfile
    try:
        conn = pyodbc.connect(conn_str)
        cursor = conn.cursor()
    except pyodbc.Error, msg:
        print("Error acessing to mdb file: %s" % msg)
        return None, None
    else:
        print("Connected succesfully to database %s" % mdbfile)
        return conn, cursor


def execute_query(cursor, QUERY, params):
    """ ищем записи """
    try:
        cursor.execute(QUERY, list(params))
    except BaseException, msg:
        print("An error ocured while executing query:", msg)
    else:
        return cursor.fetchall()


def main():
    """ ищем записи по станции в БД. stations - список кратких названий станций """
    # данные для поиска
    data = np.loadtxt(EQ_FILE, delimiter=" ", skiprows=1,
        dtype=[("num", int), ("date", "|S10"), ("time", "|S15"),
        ("lat", float), ("lon", float)])
    #
    for num, _date, _time, lat, lon in data:
        # дата/время из строк
        #dt = datetime.datetime.strptime(_date+_time, "%Y%m%d%H:%M:%S.%f")
        dt = datetime.datetime.strptime(_date+_time, "%d.%m.%Y%H:%M:%S.%f")
        # what connection to use? data before 2006 year or after?
        if dt.year < 2006: conn, cursor = conn1, cursor1
        else: conn, cursor = conn2, cursor2
        #===
        print num, _date, _time,
        sec = int( abs( (dt - DT0).total_seconds() ) ) # время (в сек) нашего события для поиска
        # params: date    time interval in seconds
        params = (dt.date(), sec-INTERVAL, sec+INTERVAL)
        # lets find idDir
        prns_dirs = execute_query(cursor, SELECT_PrnsDir, params)
        # what to do if found nothing?
        if not prns_dirs:
            print#("No events found in DB for %s" % dt)
            continue
        # always take first (0) element, but have to check it!
        idDir, Path, UTevent, DateE, TimeE, Energy = prns_dirs[0]
        # now lets find in table Prns
        params = (idDir,)
        prns = execute_query(cursor, SELECT_Prns, params)
        # now, for every prns, do job...
        timeP, timeS = None, None
        for id_prn, id_dir, fname, seisfile, K, seismgrStName in prns:
            if not seismgrStName:
                #print("No seismgrStName, but", id_prn, id_dir, fname, seisfile, K, seismgrStName)
                continue #not seisfile or seismgrStName
            for STATION in STATIONS:
                if STATION in seismgrStName.upper():
                    # волны; записано станцией station
                    waves = execute_query(cursor, SELECT_WAVES, [id_prn])
                    # output
                    for wave in waves:
                        NameWave, TimeWave, komp, Z = wave
                        if NameWave.upper().startswith("E"): continue
                        #print " ".join(map(str, wave)),
                        #
                        if NameWave.upper().startswith("P"):
                            timeP = TimeWave
                        elif NameWave.upper().startswith("S"):
                            timeS = TimeWave
                        else:
                            pass
                    print timeP, timeS,
            #==
            
        #
        print


if __name__ == "__main__":
    #===
    conn1, cursor1 = setup_db_conn(DB_FILE1)
    conn2, cursor2 = setup_db_conn(DB_FILE2)
    if conn1 is None:
        print("Error connecting to database %s" % DB_FILE1)
        sys.exit(1)
    if conn2 is None:
        print("Error connecting to database %s" % DB_FILE2)
        sys.exit(1)
    try:
        # установим соединение с бд access
        # main def
        main()
    #except BaseException, msg:
    #   print("An error ocured. Error string is:", msg)
    # закрыть соединение с бд
    finally:
        conn1.close()
        conn2.close()
    
