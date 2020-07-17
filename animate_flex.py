#! /usr/bin/env python3
# coding=utf-8
""""""
"""
Author: radenz@tropos.de


""" 

import datetime
import argparse
import sys, os
import gc
import subprocess
import numpy as np
import toml
sys.path.append("..")
import trace_source


parser = argparse.ArgumentParser()
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo')
parser.add_argument('--datetime', help='date in the format YYYYMMDD-HH')
parser.add_argument('--levels', nargs='+', type=int)
#parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')


args = parser.parse_args()


config_file = 'config_{}.toml'.format(args.station)


with open(config_file) as f:
    config = toml.loads(f.read())

end = datetime.datetime.strptime(args.datetime, '%Y%m%d-%H')
savepath = '{}/{}_maps'.format(config['plot_dir'], end.strftime('%Y%m%d_%H'))
print("savepath ", savepath)

folder = config['partposit_dir'] + '{}/'.format(end.strftime('%Y%m%d_%H'))
print('partposit_dir', folder)

dt_range = [end-datetime.timedelta(days=10), end]

files = os.listdir(folder)
files = sorted([f for f in files if 'partposit' in f])
ls = trace_source.land_sfc.land_sfc()

if args.levels is not None:
    levels = args.levels
else:
    # get levels from config file
    raise ValueError
print('levels ', args.levels)

level_to_heights = {}

for f in files[:]:

    for i in levels:
        dt = datetime.datetime.strptime(f[10:], '%Y%m%d%H%M%S')
        part_pos = trace_source.flexpart.read_partpositions(folder + f, 1, ctable=False)

        traj = trace_source.flexpart.read_flexpart_traj_meta(folder + "trajectories.txt")
        level_to_heights[i] = np.mean(traj['releases_meta'][i-1]['heights'])
        trace_source.flexpart.plot_part_loc_map(part_pos, i, dt, traj, savepath, ls=ls, config=config)
    gc.collect()


os.chdir(savepath)
print(os.getcwd())

for i in levels:

    fname_animation = "{}_{:.0f}_r{:0>2}_{}.gif".format(end.strftime('%Y%m%d_%H'), level_to_heights[i], i, args.station)
    command = "convert -scale 70% -coalesce -layers Optimize -delay 20 -loop 0 `ls r{:0>2}*.png | sort -r` {}".format(i, fname_animation)
    print('run: ', command)
    process = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)

    # from flexpart module
    # convert -scale 70% -coalesce -layers Optimize -delay 20 -loop 0 `ls r11*.png | sort -r` r11.gif
    # from notebook 
    # convert -delay 20 -loop 0 `ls r2*.png | sort -r` r2.gif
    # convert -resize 1500x1000 -delay 20 -loop 0 `ls r4*.png | sort -r` r4.gif
    # convert -scale 50% -coalesce -layers Optimize -delay 20 -loop 0 `ls r11*.png | sort -r` r11.gif
    # convert -scale 70% -coalesce -layers Optimize -delay 20 -loop 0 `ls r11*.png | sort -r` r11.gif