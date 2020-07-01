#! /usr/bin/env python3
# coding=utf-8
""""""
"""
Author: radenz@tropos.de

provide the fuctions that are needed for the generation of the input files

"""

import sys, os

import datetime
from collections import defaultdict, Counter, namedtuple
import numpy as np
import toml
#import MySQLdb
import csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../')
import trace_source


def write_hysplit_table(fname, time_list, height_list, station, runtime, coord_list=False):
    """ write a hysplit input table
    :param fname target filename
    :param station dictionary with the station information
    :param runtime backward runtime
    """

    with open(fname,'w') as f:
        f.write("YYYY MM DD HH LAT LON H_AG RUNTIME\n")

        for dt in time_list:
            if coord_list:
                print(dt, coord_list[dt.strftime("%Y%m%d-%H%M")])
                coords = coord_list[dt.strftime("%Y%m%d-%H%M")]
            else:
                coords = station['lat'], station['lon']
            for height in height_list:
                # 2017 02 26 00 34.67700 33.038 1000 -24
                f.write(dt.strftime("%Y %m %d %H") + \
                    ' {:} {:} {:>5} {:}\n'.format(*coords, height, runtime))


class gen_input():
    """
    put multiple hysplit trajectories and the statistics together in a netcdf file
    """
    def __init__(self, config_file='../config.toml'):
        """
        :param config_file toml-file with the configuration
        """

        with open(config_file) as config_file:
            self.config = toml.loads(config_file.read())

        self.config['time']['begin_dt'] = datetime.datetime.strptime(self.config['time']['begin'],
                                                                     '%Y-%m-%d_%H')
        self.config['time']['end_dt'] = datetime.datetime.strptime(self.config['time']['end'],
                                                                   '%Y-%m-%d_%H')
        print(self.config)

        self.dt_list = trace_source.time_list(self.config['time']['begin_dt'],
                                              self.config['time']['end_dt'],
                                              self.config['time']['step'])
        print('dt_list', self.dt_list)
        self.height_list = list(range(500, self.config['height']['top']+1, 500))


    def input_hysplit(self):
        """
        generate the Hysplit input table
        """

        outfile = self.config['output_dir'] \
                  + 'hysplit_input_{}_{}.csv'.format(self.config['station']['short_name'], self.dt_list[0].strftime('%Y%m%d_%H')) 

        if "moving" in self.config['station'].keys() and self.config['station']['moving'] == True:
            print("moving platform")
            track_data = {}
            with open(self.config['station']["trackfile"]) as f:
                reader = csv.reader(f, delimiter=' ')
                for row in reader:
                    dt = datetime.datetime.strptime(row[0]+" "+row[1], "%d.%m.%Y %H:%M:%S")
                    #dt = datetime.datetime.strptime(row[0]+" "+row[1], "%Y-%m-%d %H:%M:%S")
                    coords = (float(row[2]), float(row[3]))
                    track_data[dt.strftime("%Y%m%d-%H%M")] = coords

            write_hysplit_table(outfile, self.dt_list, self.height_list,
                    self.config['station'],
                    self.config['time']['tr_duration'],
                    coord_list=track_data)
            
        else:
            write_hysplit_table(outfile, self.dt_list, self.height_list,
                                self.config['station'],
                                self.config['time']['tr_duration'])
        print('hysplit input file at', outfile)


    def input_flexpart(self):
        """

        :return:
        """

        flex_svr = self.config['flexpart-server']
        flex_db = MySQLdb.connect(host=flex_svr['host'], user=flex_svr['user'],
                                  db=flex_svr['db'], passwd=flex_svr['passwd'])

        flex_cursor = flex_db.cursor()
        flex_cursor.execute("""SELECT SimulationID from TRACE_releases""")
        simids = flex_cursor.fetchall()
        print(simids)
        simids = [str(elem[0]) for elem in simids]
        print(simids)

        delta_h = 200
        delta_dt = datetime.timedelta(minutes=60)
        sim_dur = datetime.timedelta(hours= abs(self.config['time']['tr_duration']))

        lat = self.config['station']['lat']
        lon = self.config['station']['lon']
        comment = 'trace_source cyprus'

        i = 0
        for dt in self.dt_list:
            for ih, height in enumerate(self.height_list):
                simID = dt.strftime('%y%m%d') + '{:0>3}'.format(i)
                bottom = height - delta_h
                top = height + delta_h
                startsim = (dt-sim_dur).strftime('%Y-%m-%d %H:%M:%S')
                startt = (dt-delta_dt).strftime('%Y-%m-%d %H:%M:%S')
                stopt = (dt+delta_dt).strftime('%Y-%m-%d %H:%M:%S')
                stopsim = stopt

                if not simID in simids:
                    # add the releases
                    flex_cursor.execute("SELECT SimulationID from TRACE_releases")
                    flex = flex_cursor.fetchall()
                    sqlcmd = """INSERT INTO TRACE_releases
                    (ID, SimulationID, comment, Latitude, Longitude, Bottom, Top, startRelease, StopRelease)
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"""
                    print((str(len(flex) + 1), simID, comment, lat, lon, bottom, top, startt, stopt))
                    flex_cursor.execute(sqlcmd,
                                        (str(len(flex) + 1), simID, comment, lat, lon,
                                         bottom, top, startt, stopt))

                    # add the simulations
                    flex_cursor.execute("SELECT ID from TRACE_simulations")
                    flex = flex_cursor.fetchall()
                    sqlcmd = """INSERT INTO TRACE_simulations
                    (ID, MeasurementName, StartSim, StopSim, scheduled, Server)
                    VALUES (%s, %s, %s, %s, %s, %s)"""
                    print((simID, comment, startsim, stopsim, 127, flex_svr['host']))
                    flex_cursor.execute(sqlcmd, (simID, comment, startsim, stopsim, 127, flex_svr['host']))

                else:
                    print('!! simid already there')
                i += 1


        flex_cursor.close()
        flex_db.rollback()
        flex_db.close()


if __name__ == '__main__':
    #tr_input = gen_input(config_file='config.toml')
    #tr_input.input_flexpart()

    # tr_input = gen_input(config_file='config_barbados.toml')
    # tr_input.input_hysplit()

    tr_input = gen_input(config_file='config.toml')
    #tr_input = gen_input(config_file='config_ps106.toml')
    tr_input.input_hysplit()
