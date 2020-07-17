#!/usr/bin/env python3

#
# radenz@tropos.de
# 2019-09-23 inital file on lidardaten
# 2020-07-07 modified to work with the docker setup
#

import os, shutil
import datetime
import argparse
import subprocess
import string
import toml

def gen_dt_list(begin, end, step):
    print('daterange ', begin, 'to', end)
    dt_list = []
    current = begin
    while current <= end:
        dt_list.append(current)
        current += datetime.timedelta(hours=step)
    return dt_list

def refresh_AVAILABLE(data_dir):

    files = os.listdir(data_dir)
    files.sort()
    with open('AVAILABLE', 'w') as f:
        f.write('DATE     TIME         FILENAME     SPECIFICATIONS' + '\n')
        f.write('YYYYMMDD HHMISS' + '\n')
        f.write('________ ______      __________      __________' + '\n')
        for file in files:
            if file.isdigit() :#<> 'A') and (file[0] <> '-'):
                date = file[0:8]+ ' '
                time = file[8:10] + '0000      '
                f.write(date + time + file + '      ' + 'ON_DISC' + '\n')
    print("updated AVAILABLE")



def write_pathnames(control_dir, out_dir, data_dir, available_file):
    
    with open('pathnames_template') as t_file:
        temp = t_file.read()
      
    t = string.Template(temp)
    # all pathes with tailing /
    pathnames = t.substitute(control_dir=control_dir,
                             out_dir=out_dir,
                             data_dir=data_dir,
                             available_file=available_file)
    #print(pathnames)
    with open('pathnames', 'w') as file:
        file.write(pathnames)
    print("written pathnames")



def write_COMMAND(control_dir, begin, end):
    
    with open('COMMAND_template') as t_file:
        temp = t_file.read()
      
    t = string.Template(temp)
    out_string = t.substitute(begin_date=begin.strftime("%Y%m%d"),
                              begin_time=begin.strftime("%H%M%S"),
                              end_date=end.strftime("%Y%m%d"),
                              end_time=end.strftime("%H%M%S"))
    #print(pathnames)
    with open(control_dir + 'COMMAND', 'w') as file:
        file.write(out_string)
    print("written COMMAND")

    

def write_RELEASES(control_dir, time, bbox, heights, no_part, comment):
    """

    """

    
    with open('RELEASES_header_template') as t_file:
        out_string = t_file.read()

    with open('RELEASES_template') as t_file:
        temp = t_file.read()
    t = string.Template(temp)

    for h in heights:
        c = comment + str((h[0]+h[1])/2.) + 'm ' + str(no_part)
        out_string += t.substitute(begin=time[0].strftime("%Y%m%d %H%M%S"),
                                   end=time[1].strftime("%Y%m%d %H%M%S"),
                                   lon_ll = '%9.4f' % bbox['lon_ll'],
                                   lat_ll = '%9.4f' % bbox['lat_ll'],
                                   lon_ur = '%9.4f' % bbox['lon_ur'],
                                   lat_ur = '%9.4f' % bbox['lat_ur'],
                                   lower = '%9.3f' % h[0],
                                   upper = '%9.3f' % h[1],
                                   no_part = '%9d' % no_part,
                                   comment = c
                                  )
    #print(pathnames)
    with open(control_dir + 'RELEASES', 'w') as file:
        file.write(out_string)
    print("written RELEASES")





parser = argparse.ArgumentParser(usage="usage: %prog [date-interval]")
#parser.add_option("interval", help="interval %Y%m%d-%Y%m%d")
parser.add_argument("-f", "--fields",
                  default='gfs1',
                  help="select the model data {gfs1, gfs0p25}")
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo')
#parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')
parser.add_argument('--date', help='date in the format YYYYMMDD')
args = parser.parse_args()

print('args', args)

day = datetime.datetime.strptime(args.date, "%Y%m%d")

dt_list = gen_dt_list(day, day + datetime.timedelta(hours=21), 3)

print(dt_list)


config = toml.load('config_{}.toml'.format(args.station))
# station = 'punta'
# coord = [-53.1344, -70.8801]
# #station = 'limassol'
# #coord = [34.677, 33.038]

# station = 'PS113'
# coord = [19.8, -21.2]

# #station = 'PS122'
# #coord = [85.1, 133.7]
# #coord = [84.9, 136.0]
# #coord = [86.1, 113.1]
# #coord = [84.8, 134.8] #20191013

comment = 'trace {} '.format(args.station)

if args.fields == 'gfs1':
    data_dir = "../data/gfs_083.2/"
elif args.fields == 'gfs0p25':
    data_dir = "../data/gfs_083.2/"
else:
    raise ValueError

out_dir = '../flexpart_partposit/{}/'.format(args.station)
print('out_dir', out_dir)

os.chdir('flexpart_simulations')

if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

refresh_AVAILABLE(data_dir)

write_pathnames(os.getcwd(),
                out_dir+"out/",
                data_dir,
                "AVAILABLE")

for dt in dt_list[:]:
    print(dt)
    if not os.path.isdir(out_dir+"out/"):
        os.mkdir(out_dir+"out/")
    write_COMMAND(os.getcwd()+'/', dt-datetime.timedelta(days=10), dt)

    time = [dt-datetime.timedelta(minutes=5), dt]
    bbox = {'lon_ll': config['station']['lon']-0.1, 'lat_ll': config['station']['lat']-0.1,
            'lon_ur': config['station']['lon']+0.1, 'lat_ur': config['station']['lat']+0.1}

    center_heights = list(range(500, config['height']['top']+1, 500))
    heights = [(h-200, h+200) for h in center_heights]
    #assert heights == [(300, 700), (800,1200), (1300, 1700), (1800,2200),
    #           (2300, 2700), (2800,3200), 
    #           (3300, 3700), (3800,4200), (4300, 4700), (4800,5200),
    #           (5300, 5700), (5800,6200), (6300, 6700), (6800,7200),
    #           (7300, 7700), (7800,8200), (8300, 8700), (8800,9200),
    #           (9300, 9700), (9800,10200), (10300, 10700), (10800,11200),
    #           (11300, 11700), (11800,12200), 
    #           ]
    #heights = [(300,700)]
    write_RELEASES(os.getcwd()+'/', time, bbox, heights, config['flexpart']['no_particles'], comment)

    # if options.model == 'gfs1':
    #     flexpart_bin = ['/home/flexpart/flexpart/programs/flexpart/FLEXPART_GFS', '.']
    #     # thesting the flexpart version 9.2 as suggested in the issue
    #     #flexpart_bin = ['/home/flexpart/flexpart/programs/flexpart/FLEXPART_GFS_highres', '-v2']
    # elif options.model == 'gfs0p25':
    #     flexpart_bin = ['/home/flexpart/flexpart/programs/flexpart/FLEXPART_GFS_highres', '-v']
    # else:
    #     raise ValueError

    #p = subprocess.Popen(flexpart_bin, bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print('cwd ', os.getcwd())
    #print('content dir ', os.listdir('.'))
    process = subprocess.run('FLEXPART', shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True)
    print(process.stdout)

    print('remove and rename ', out_dir+"out", out_dir+dt.strftime("%Y%m%d_%H"))
    if dt.strftime("%Y%m%d_%H") in os.listdir(out_dir):
        shutil.rmtree(out_dir+dt.strftime("%Y%m%d_%H"))
    os.rename(out_dir+"out", out_dir+dt.strftime("%Y%m%d_%H"))

