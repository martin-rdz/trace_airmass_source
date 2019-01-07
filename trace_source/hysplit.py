#! /usr/bin/env python3
# coding=utf-8
""""""
"""
Author: radenz@tropos.de

TODO
----

""" 

import sys, os
import re
import gc
import datetime
from collections import defaultdict, Counter, namedtuple
import numpy as np
import toml
import netCDF4

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + '/../')
import trace_source


def plot_trajectories_ens(traj, savepath, ls=None, config=None):
    """
    plot multiple trajectories into one scene

    Args:
        traj (:class:`.trajectory`): trajectory to plot instance to plot
        savepath (str): path to save
        ls (:class:`trace_source.land_sfc.land_sfc`, optional): pre loaded land surface information
        config (dict, optional): the toml derived config dict

    Returns:
        None
    """

    import matplotlib
    matplotlib.use('Agg')
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    if ls is None:
        ls = trace_source.land_sfc.land_sfc()

    colors = ['purple', 'darkorange', 'hotpink', 'dimgrey']
    c = 'purple'

    if not os.path.isdir(savepath):
        os.makedirs(savepath)

    ###
    # the map
    ###

    fig = plt.figure(figsize=(8, 10))
    if config is not None and "centerlon" in config['plotmap']:
        ax = plt.axes(projection=ccrs.Miller(central_longitude=config['plotmap']['centerlon']))
    else:
        ax = plt.axes(projection=ccrs.Miller(central_longitude=-170.))
        ax = plt.axes(projection=ccrs.Miller())
        # Projection for the north pole
        # ax = plt.axes(projection=ccrs.NorthPolarStereo())
        raise ValueError('provide plotmap.centerlon in the config file')

    ####
    # make a color map of fixed colors
    cmap = matplotlib.colors.ListedColormap(['lightskyblue', 'darkgreen', 'khaki', 'palegreen', 'red', 'white', 'tan'])
    bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    ####
    pcm = ax.pcolormesh(ls.longs, ls.lats, ls.land_sfc, cmap=cmap, norm=norm, transform=ccrs.PlateCarree())
    ax.coastlines()

    for k,v in traj.data.items():
        ax.plot(v['longitude'], v['latitude'],
                linewidth=1.5,
                color=c,
                transform=ccrs.Geodetic())

        ax.plot(v['longitude'][::24], v['latitude'][::24], '.',
                color=c,
                transform=ccrs.Geodetic())

    ax.gridlines(linestyle=':')
    if config is not None and "bounds" in config['plotmap']:
        ax.set_extent(config['plotmap']['bounds'], crs=ccrs.PlateCarree())
    else:
        # bounds for punta arenas
        ax.set_extent([-180, 180, -75, 0], crs=ccrs.PlateCarree())
        # bounds for atlantic
        #ax.set_extent([-180, 80, -70, 75], crs=ccrs.PlateCarree())
        # bounds for the north pole
        #ax.set_extent([-180, 180, 50, 90], crs=ccrs.PlateCarree())
        raise ValueError('provide plotmap.bounds in the config file')
    ax.set_title('Trajectory {}UTC\n{} {} '.format(traj.info['date'].strftime('%Y-%m-%d %H'),
                                                   traj.info['lat'], traj.info['lon']),
                 fontweight='semibold', fontsize=13)

    savename = savepath + "/" + traj.info['date'].strftime("%Y%m%d_%H") \
        + '_' + '{:0>5.0f}'.format(traj.info['height']) + "_trajectories_map.png"
    fig.savefig(savename, dpi=250)
    plt.close()

    ###
    # the profile
    ###

    fig, ax = plt.subplots(4, 1, figsize=(6, 7), sharex=True)

    ax[0].grid(True, axis='y')
    for k,v in traj.data.items():
        ax[0].plot(v['time'], v['height'])
        ax[1].plot(v['time'], v['AIR_TEMP'][:] - 273.15)
        ax[2].plot(v['time'], v['RELHUMID'][:])
        ax[3].plot(v['time'], v['RAINFALL'][:])

    ax[0].set_ylim(0, 10000)
    ax[0].set_ylabel('Height [m]')
    ax[0].xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=2))
    ax[0].xaxis.set_minor_locator(matplotlib.dates.HourLocator([0, 12]))
    ax[0].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))

    # ax[1].set_ylim(0,10000)
    ax[1].set_ylabel('Temperature [°C]')
    ax[1].set_ylim(-40, 10)
    ax[1].xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=2))
    ax[1].xaxis.set_minor_locator(matplotlib.dates.HourLocator([0, 12]))
    ax[1].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))

    # ax[2].grid(True, axis='y')
    ax[2].set_ylim(0, 100)
    ax[2].set_ylabel('rel Humidity [%]')
    ax[2].xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=2))
    ax[2].xaxis.set_minor_locator(matplotlib.dates.HourLocator([0, 12]))
    ax[2].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))

    # ax[3].grid(True, axis='y')
    ax[3].set_ylim(0, 20)
    ax[3].set_xlabel('Day')
    ax[3].set_ylabel('Precip [mm]')
    ax[3].xaxis.set_major_locator(matplotlib.dates.DayLocator(interval=2))
    ax[3].xaxis.set_minor_locator(matplotlib.dates.HourLocator([0, 12]))
    ax[3].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%m-%d'))

    for axis in ax:
        axis.tick_params(axis='both', which='both', right=True, top=True,
                         direction='in')

    ax[0].set_title('Trajectory {}UTC\n{} {} '.format(traj.info['date'].strftime('%Y-%m-%d %H'),
                                                      traj.info['lat'], traj.info['lon']),
                    fontweight='semibold', fontsize=13)
    fig.tight_layout()

    savename = savepath + "/" + traj.info['date'].strftime("%Y%m%d_%H") \
        + '_' + '{:0>5.0f}'.format(traj.info['height']) + "_trajectories_prof.png"
    fig.savefig(savename, dpi=250)
    plt.close()


def save_item(dataset, item_data):
    """
    Save an item to the dataset with the data given as a dict

    Args:
        dataset (:obj:netCDF4.Dataset): netcdf4 Dataset to add
        item_data (dict): with the data to add, for example:

    ==================  ===============================================================
     Key                 Example                            
    ==================  ===============================================================
     ``var_name``        Z                                  
     ``dimension``       ('time', 'height')                 
     ``arr``             self.corr_refl_reg[:].filled()     
     ``long_name``       "Reflectivity factor"              
     **optional**                                             
     ``comment``         "Wind profiler reflectivity factor corrected by cloud radar"
     ``units``           "dBz"                              
     ``missing_value``   -200.                              
     ``plot_range``      [-50., 20.]                        
     ``plot_scale``      "linear"                           
     ``vartype``         np.float32                         
    ==================  ===============================================================
 
    """

    if 'vartype' in item_data.keys():
        item = dataset.createVariable(item_data['var_name'], item_data['vartype'], item_data['dimension'])
    else:
        item = dataset.createVariable(item_data['var_name'], np.float32, item_data['dimension'])
    item[:] = item_data['arr']
    item.long_name = item_data['long_name']
    if 'comment' in item_data.keys():
        item.comment = item_data['comment']
    if 'units' in item_data.keys():
        item.units = item_data['units']
    if 'units_html' in item_data.keys():
        item.units_html = item_data['units_html']
    if 'missing_value' in item_data.keys():
        item.missing_value = item_data['missing_value']
    if 'plot_range' in item_data.keys():
        item.plot_range = item_data['plot_range']
    if 'plot_scale' in item_data.keys():
        item.plot_scale = item_data['plot_scale']
    if 'axis' in item_data.keys():
        item.axis = item_data['axis']

    return dataset


class assemble_time_height():
    """
    assemble a time height period by putting multiple hysplit trajectories togehter,
    calculate the statistics and save to a netcdf file
    
    Args:

        config_file (str, optional): path to the config file 

    """
    def __init__(self, config_file='../config.toml'):
        with open(config_file) as config_file:
            self.config = toml.loads(config_file.read())

        self.config['time']['begin_dt'] = datetime.datetime.strptime(self.config['time']['begin'],
                                                                     '%Y-%m-%d_%H')
        self.config['time']['end_dt'] = datetime.datetime.strptime(self.config['time']['end'],
                                                                   '%Y-%m-%d_%H')
        print('config', self.config)
        self.dt_list = trace_source.time_list(self.config['time']['begin_dt'],
                                              self.config['time']['end_dt'],
                                              self.config['time']['step'])
        print('dt_list', self.dt_list)
        self.height_list = list(range(500, self.config['height']['top']+1, 500))


    def assemble(self, dt_range=None):
        """
        assemble the statistics for a range of trajectories and
        save the statistics to dicts
        
        Args:
            dt_range (list(datetime), optional): timerange for that the statistics is assembled,
                default taken from config 

        """
        if dt_range is not None:
            self.dt_list = trace_source.time_list(dt_range[0],
                                                  dt_range[1],
                                                  self.config['time']['step'])

        # only for the testcase
        traj_dir = self.config['traj_dir']
        files = os.listdir(traj_dir)
        # filter only for the trajectory files with tdump extension
        files = [f for f in files if f[-6:] == '.tdump']

        # the defaultdict is used here to sort the files by datetime within a dictionary
        filtered_files = defaultdict(list)
        for f in files:
            # regex the yyyymmdd-hh timestamp in the filename
            dt = datetime.datetime.strptime(re.search('([0-9]{8})-([0-9]){2}', f).group(0), '%Y%m%d-%H')
            height = float(re.search('([0-9]{3,6})(?=_0[0-9-]{1,4}.tdump)', f).group(0))
            #print(f, dt, height)
            if dt >= self.dt_list[0] and dt <= self.dt_list[-1]:
                filtered_files[dt].append((f,height))

        # here an empty dict is generated with a zero containing array
        self.stat2d_dict = defaultdict(lambda: np.zeros((len(self.dt_list), len(self.height_list))))
        self.statls_dict = defaultdict(lambda: np.zeros((len(self.dt_list), len(self.height_list), 7)))

        self.raw_dict = defaultdict(lambda: np.zeros((len(self.dt_list), len(self.height_list),
                                                     abs(self.config['time']['tr_duration'])+1)))

        # TODO make more than 7 geo names possible
        ng = trace_source.land_sfc.named_geography(self.config['geonames'])
        self.geo_names = ng.geo_names
        no_geo_names = len(list(self.geo_names.keys()))
        self.statgn_dict = defaultdict(lambda: np.zeros((len(self.dt_list),
                                                         len(self.height_list),
                                                         no_geo_names)))

        ls = trace_source.land_sfc.land_sfc()
        self.ls_categories = ls.categories

        for it, dt in enumerate(self.dt_list[:]):
            print(dt)
            # sort by height
            f_list = sorted(filtered_files[dt], key= lambda x: x[1])
            print('file_list ', f_list)
            assert len(f_list) > 1
            for ih, f in enumerate(f_list):
                print(it, ih, f[1], dt)
                traj = trajectory(self.config)
                traj.load_file(traj_dir+f[0], silent=True)
                savepath = '{}/{}'.format(self.config['plot_dir'], dt.strftime('%Y%m%d'))


                if "timeinterval" in self.config['plotmap']:
                    timeinterval = self.config['plotmap']['timeinterval']
                else:
                    timeinterval = 12
                if "heights" in self.config['plotmap']:
                    heightlist = self.config['plotmap']['heights']
                else:
                    heightlist = [1500.0, 3000.0, 4500.0]
                #if f[1] == 3000.0 and dt.hour % 12 == 0:
                if f[1] in heightlist and dt.hour % timeinterval == 0:
                    print("plotting ", f[1], dt.hour)
                    plot_trajectories_ens(traj, savepath, ls=ls, config=self.config)
                #continue

                traj.evaluate(silent=True)
                traj.add_land_sfc(ls, silent=True)
                traj.add_ensemble_land_sfc(ls)
                traj.add_ensemble_geo_names(ng)
                #traj.add_area_land_sfc('md', ls, silent=True)
                #traj.add_area_land_sfc(2000, ls, silent=True)

                #print("at step", it, dt, ih, f)
                #print('keys ', traj.statistics.keys())
                # now the empty dict is filled with the keys (and values) of the statistics dict from traj
                for k in list(traj.statistics.keys()):
                    self.stat2d_dict[k][it, ih] = traj.statistics[k]
                # subset of trajectory data to collect
                param_collect = ['latitude', 'longitude', 'height', "PRESSURE", "AIR_TEMP",
                                 "RAINFALL", "RELHUMID", "TERR_MSL", 'age']
                if 'land_sfc_category' in list(traj.data.keys()):
                    param_collect.append('land_sfc_category')
                for k in param_collect:
                    #self.raw_dict[k][it, ih, :traj.data[1][k].shape[0]] = traj.data[1][k]
                    self.raw_dict[k][it, ih, :] = traj.data[1][k]
                    #self.raw_dict[k][it, ih, traj.data[1][k].shape[0]:] = -999.

                for k in list(traj.stat_ls.keys()):
                    self.stat2d_dict[k+'_no_below'][it, ih] = traj.stat_ls[k].no_below
                    print('stat ls ', k, traj.stat_ls[k])
                    self.statls_dict[k][it, ih] = list(traj.stat_ls[k].counter.values())

                for k in list(traj.stat_gn.keys()):
                    self.stat2d_dict[k+'_no_below'][it, ih] = traj.stat_gn[k].no_below
                    print('stat gn ', k, traj.stat_gn[k])
                    self.statgn_dict[k][it, ih] = list(traj.stat_gn[k].counter.values())

        # trying to free memory
        del ls
        del ng


    def dump2netcdf(self):
        """
        dump the assembled data to a netcdf file repeatatly call :meth:`save_item`
        configuration (directories, names, etc) is given by the config file
        """

        timestamps = np.array([(dt - datetime.datetime(1970,1,1)).total_seconds() for dt in self.dt_list])
        hours_cn = np.array([dt.hour + dt.minute / 60. + dt.second / 3600. for dt in self.dt_list])

        if not os.path.isdir(self.config['output_dir']):
            os.makedirs(self.config['output_dir'])

        ncfile = self.config['output_dir'] +\
                '{}_{}_hysplit-output.nc'.format(self.dt_list[0].strftime('%Y%m%d'), self.config['station']['short_name'])
        # ncfile = "/home/devel/" +\
        #          '{}_hysplit_output.nc'.format(self.config['time']['begin_dt'].strftime('%Y%m%d'))

        #dataset = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
        dataset = netCDF4.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')

        dim_time = dataset.createDimension('time', len(self.dt_list))
        dim_height = dataset.createDimension('height', len(self.height_list))
        dim_age = dataset.createDimension('time_age', abs(self.config['time']['tr_duration'])+1)
        dim_cat = dataset.createDimension('categories', 7)
        dim_regions = dataset.createDimension('regions', len(list(self.geo_names.keys())))

        # times_cn = dataset.createVariable('time', np.float32, ('time',))
        # times_cn[:] = hours_cn.astype(np.float32)
        # times_cn.units = "hours since " + self.begin_dt.strftime('%Y-%m-%d') + "00:00:00 +00:00"
        # times_cn.long_name = "Decimal hours from midnight UTC"
        # times_cn.axis = "T"

        save_item(dataset, {'var_name': 'timestamp', 'dimension': ('time', ),
                            'vartype': 'i4',
                            'arr': timestamps.astype(np.int32),
                            'long_name': "Unix timestamp",
                            'units': "s", 'axis': 'T'})
        save_item(dataset, {'var_name': 'time', 'dimension': ('time', ),
                            'arr': hours_cn.astype(np.float32),
                            'long_name': "Decimal hours from midnight UTC",
                            'units': "hours since {} 00:00:00 +00:00".format(self.dt_list[0].strftime('%Y-%m-%d')),
                            'axis': 'T'})
        save_item(dataset, {'var_name': 'range', 'dimension': ('height', ),
                            'arr': np.array(self.height_list).astype(np.float32)/1000.,
                            'long_name': "Height",
                            'units': "km", 'axis': 'Z'})
        save_item(dataset, {'var_name': 'age', 'dimension': ('time_age', ),
                            'arr': self.raw_dict['age'][0, 0],
                            'long_name': "Age of trajectory",
                            'units': "h"})

        for k in list(self.stat2d_dict.keys()):
            print(k, self.stat2d_dict.get(k).shape)
            dataset = save_item(dataset, {'var_name': k, 'dimension': ('time', 'height'),
                                          'arr': self.stat2d_dict.get(k).copy().astype(np.float32),
                                          'long_name': k})

        # its sufficient to save the age of the trajectory once
        raw_data_keys = list(self.raw_dict.keys())
        raw_data_keys.remove('age')
        # chance to modify some parameter descriptions for better readability
        modified_params = {key: {'var_name': key.lower(),
                                 'long_name': "Hysplit " + key.lower()} for key in raw_data_keys}
        modified_params['height'] = {'var_name': 'traj_height', 'long_name': "Hysplit height of air parcel"}
        if 'land_sfc_category' in list(modified_params.keys()):
            modified_params['land_sfc_category']['long_name'] = "Modis land use category (simplified)"
        for k in raw_data_keys:
            print(k, self.raw_dict.get(k).shape)
            dataset = save_item(dataset, {'var_name': modified_params[k]['var_name'],
                                          'dimension': ('time', 'height', 'time_age'),
                                          'arr': self.raw_dict.get(k).copy().astype(np.float32),
                                          'long_name': modified_params[k]['long_name']})

        # save the land use
        ls_data_keys = list(self.statls_dict.keys())
        print('ls_data_keys ', ls_data_keys)
        modified_params = {key: {'var_name': key,
                                 'long_name': "land use " + key.lower(),
                                 'comment': str(self.ls_categories)} for key in ls_data_keys}
        for k in ls_data_keys:
            print(k, self.statls_dict.get(k).shape)
            dataset = save_item(dataset, {'var_name': modified_params[k]['var_name'],
                                          'dimension': ('time', 'height', 'categories'),
                                          'arr': self.statls_dict.get(k),
                                          'long_name': modified_params[k]['long_name'],
                                          'comment': modified_params[k]['comment']})


        # and the geo names
        gn_data_keys = list(self.statgn_dict.keys())
        print('gn_data_keys ', gn_data_keys)
        modified_params = {key: {'var_name': key,
                                 'long_name': "geography names " + key.lower(),
                                 'comment': str(self.geo_names)} for key in gn_data_keys}
        for k in gn_data_keys:
            print(k, self.statgn_dict.get(k).shape)
            dataset = save_item(dataset, {'var_name': modified_params[k]['var_name'],
                                          'dimension': ('time', 'height', 'regions'),
                                          'arr': self.statgn_dict.get(k),
                                          'long_name': modified_params[k]['long_name'],
                                          'comment': modified_params[k]['comment']})

        # save_item(dataset, {'var_name': 'width', 'dimension': ('time', 'height'),
        #                     'arr': .corr_width_reg[:].filled(), 'long_name': "Spectral width",
        #                     'comment': "Wind profiler spectral width (standard deviation) corrected by cloud radar (only Bragg contribution)",
        #                     'units': "m s-1", 'units_html': "m s<sup>-1</sup>",
        #                     'missing_value': -99., 'plot_range': (0.01, 4.),
        #                     'plot_scale': "logarithmic"})


        with open('output_meta.toml') as output_meta:
            meta_info = toml.loads(output_meta.read())

        dataset.description = "trace_source trajectory"
        dataset.location = self.config['station']['name']
        if "moving" in self.config['station'].keys() and self.config['station']['moving'] == True:
            dataset.coordinates = "Moving Platform!"
        else:
            dataset.coordinates = (self.config['station']['lat'], self.config['station']['lon'])
        dataset.institution = meta_info["institution"]
        dataset.authors = meta_info["authors"]
        dataset.contact = meta_info["contact"]
        dataset.creation_time = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
        dataset.day = self.dt_list[0].day
        dataset.month = self.dt_list[0].month
        dataset.year = self.dt_list[0].year
        dataset.close()
        gc.collect()


class trajectory():
    """
    handle a single trajectory
    
    Args:
        config (dict): content of the config file
    """

    def __init__(self, config):
        self.statistics = {}
        self.stat_ls = {}
        self.stat_gn = {}
        self.shapes = {}
        self.config = config
        
        
    def load_file(self, filename, silent=False):
        """load a trajectory from the .tdump file

        Args:
            filename (str): path and name of the file
            silent (bool, optional): verbose output or not
        
        """
        print('filename', filename) if not silent else None

        #parameters_key = ["PRESSURE", "AIR_TEMP", "RAINFALL", "RELHUMID", "TERR_MSL"]
        parameters_key = ["PRESSURE", "AIR_TEMP", "RAINFALL", "MIXDEPTH", "RELHUMID", "TERR_MSL"]

        # empty dict for the ensemble trajectories
        data = defaultdict(lambda: {'time': [], 'latitude': [], 'longitude': [],
                                    'age': [], 'height': [], 'parameters': []})
        
        with open(filename) as f:
            header = []
            while True:
                line = f.readline().strip()
                header.append(line)
                if 'PRESSURE' in line:
                    break
            #header = [f.readline().strip() for i in range(7)]
            int_info = int(header[0].split()[0]) + 2
            no_traj = int(header[int_info-1].split()[0])

            arr = np.empty((no_traj, 241))
            arr.fill(-999.)
            arrp = np.empty((no_traj, 241, 6))
            arrp.fill(-999.)
            arr_data = {'latitude':arr.copy(), 'longitude':arr.copy(),
                        'height':arr.copy(), 'parameters':arrp.copy()}

            idx_traj = 1
            for row in f:
                row_contents = row.split()
                # print(row_contents)
                k = int(row_contents[0])
                t = int(abs(float(row_contents[8])))
                # fill missing steps

                #print(datetime.datetime.strptime('-'.join(row_contents[2:7]), '%y-%m-%d-%H-%M'))
                #data[k]['time'].append(datetime.datetime.strptime('-'.join(row_contents[2:7]), '%y-%m-%d-%H-%M'))
                #print('age', row_contents[8])
                arr_data['latitude'][k-1, t] = float(row_contents[9])
                arr_data['longitude'][k-1, t] = float(row_contents[10])
                arr_data['height'][k-1, t] = float(row_contents[11])
                #print([float(elem) for elem in row_contents[12:]])
                arr_data['parameters'][k-1, t, :] = [float(elem) for elem in row_contents[12:]]


        infostring = header[int_info].split()
        year = filename.split('-')[1][0:4]
        self.info = {'date': datetime.datetime(int(year), int(infostring[1]),
                                               int(infostring[2]), int(infostring[3])),
                     'lat': float(infostring[4]), 'lon': float(infostring[5]),
                     'height': float(infostring[6]), 'no_traj': no_traj}

        age = np.arange(0,-241,-1)
        time = [self.info['date']+datetime.timedelta(hours=h) for h in age.tolist()]
        self.data = {}
        for k in range(1,28):
            self.data[k] = {'time': time, 'age': age,
                            'latitude': np.ma.masked_equal(arr_data['latitude'][k-1, :], -999.),
                            'longitude': np.ma.masked_equal(arr_data['longitude'][k-1, :], -999.),
                            'height': np.ma.masked_equal(arr_data['height'][k-1, :], -999.)}

            for i, param in enumerate(parameters_key):
                self.data[k][param] = np.ma.masked_equal(arr_data['parameters'][k-1, :, i], -999.)

        print('trajectory data') if not silent else None
        for k,v in self.data.items():
            print(k, list(v.keys())) if not silent else None


    def load_dict(self, d):
        """load the trajectory without the original data
        
        .. note::
            only the first member of a ensemble is saved

        Args:
            d (dict)
        """

        self.info = d['info']
        self.data = {1: d['data']}


    def load_netcdf(self, ncf, dt, height):
        """
        load the trajectory at a given datetime and height from the provided netcdf file

        Args:
            ncf (:obj:netCDF4.Dataset): netcdf4 Dataset to add
            dt: datetime of trajectory to load
            height: height of trajectory to load

        Returns:
            dict with data

            Keys: ``time``, ``age``, ``latitude``, ``longitude``, ``height``

            Optional ``pressure``, ``air_tem``, ``rainfall``, ``relhumid``, ``terr_msl``
        """

        d = {}

        d['info'] = {'date': dt, 'height': height,
                     'lat': ncf.coordinates[0], 'lon': ncf.coordinates[1]}

        dt_list = [datetime.datetime.fromtimestamp(time) for time in ncf.variables["timestamp"][:]]
        it = dt_list.index(dt)
        ih = np.where(ncf.variables["range"][:] == height)[0]

        traj_time = [dt - datetime.timedelta(hours=abs(h)) for h in ncf.variables['age'][:].tolist()]
        d['data'] = {'time': traj_time, 'age': ncf.variables['age'][:],
                     'latitude': ncf.variables['latitude'][it, ih, :][0],
                     'longitude': ncf.variables['longitude'][it, ih, :][0],
                     'height': ncf.variables['traj_height'][it, ih, :][0]}

        parameters_key = ["PRESSURE", "AIR_TEMP", "RAINFALL", "RELHUMID", "TERR_MSL"]
        for i, param in enumerate(parameters_key):
            d['data'][param] = ncf.variables[param.lower()][it, ih, :][0]

        self.load_dict(d)
        return d


    def add_land_sfc(self, ls=None, silent=False):
        """ add the land surface information to a single trajectory
            
        Args:
            ls (:class:`trace_source.land_sfc.land_sfc`, optional): pre loaded land surface information
                (separate loading consumes lot of time)
            silent (bool, optional): verbose output or not

        """

        if ls is None:
            ls = trace_source.land_sfc.land_sfc()


        occ_stat = namedtuple('occ_stat', 'no_below counter')
        for rh in self.config['height']['reception']:
            category = ls.get_land_sfc(self.data[1]['latitude'], self.data[1]['longitude'])
            print(category) if not silent else None

            if rh == 'md':
                cat = category[self.data[1]['height'] < self.data[1]['MIXDEPTH']]
            else:
                cat = category[self.data[1]['height'] < float(rh)*1000]

            no = float(cat.shape[0]) if cat.shape[0] > 0 else -1
            c = {x: cat.tolist().count(x)/float(no) for x in list(ls.categories.keys())}
            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh
            self.stat_ls['occ_below'+rh_string] = occ_stat(no_below=no, counter=c)

        print('occ_one ', occ_belowmd) if not silent else None
        print('occ_one ', occ_below2km) if not silent else None


    def add_ensemble_land_sfc(self, ls=None):
        """ add the land surface information to an ensemble trajectory
            
        Args:
            ls (:class:`trace_source.land_sfc.land_sfc`, optional): pre loaded land surface information
                (separate loading consumes lot of time)
            silent (bool, optional): verbose output or not

        """

        if ls is None:
            ls = trace_source.land_sfc.land_sfc()

        if self.info['no_traj'] == 1:
            raise ValueError('tried to call add_ensemble_land_sfc with single traj')

        occ_stat = namedtuple('occ_stat', 'no_below counter')
        for rh in self.config['height']['reception']:
            categories = np.empty((0,))
            for k, v in self.data.items():
                category = ls.get_land_sfc(v['latitude'],
                                           v['longitude'])
                if rh == 'md':
                    cat = category[v['height'] < v['MIXDEPTH']]
                else:
                    cat = category[v['height'] < float(rh)*1000]
                categories = np.append(categories, cat)

            no = float(categories.shape[0]) if categories.shape[0] > 0 else -1
            c = {x: categories.tolist().count(x)/float(no) for x in list(ls.categories.keys())}
            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh
            self.stat_ls['occ_ens_below' + rh_string] = occ_stat(no_below=no, counter=c)


    def add_ensemble_geo_names(self, ng=None):
        """
        add the geographical names to a ensemble trajectory
            
        Args:
            ng (:class:`trace_source.land_sfc.named_geography`, optional): pre loaded named geography information
                (separate loading consumes lot of time)
            silent (bool, optional): verbose output or not

        """

        if ng is None:
            ng = trace_source.land_sfc.named_geography(self.config['geonames'])

        if self.info['no_traj'] == 1:
            raise ValueError('tried to call add_ensemble_land_sfc with single traj')

        occ_stat = namedtuple('occ_stat', 'no_below counter')
        for rh in self.config['height']['reception']:
            categories = np.empty((0,))
            for k, v in self.data.items():
                category = ng.get_geo_names(v['latitude'], 
                                        v['longitude']) 
                if rh == 'md':
                    cat = category[v['height'] < v['MIXDEPTH']]
                else:
                    cat = category[v['height'] < float(rh)*1000]
                categories = np.append(categories, cat)

            no = float(categories.shape[0]) if categories.shape[0] > 0 else -1
            c = {x: categories.tolist().count(x)/float(no) for x in list(ng.geo_names.keys())}
            if rh != 'md':
                rh_string = rh + 'km'
            else:
                rh_string = rh
            self.stat_gn['region_ens_below' + rh_string] = occ_stat(no_below=no, counter=c)


    def add_area_land_sfc(self, maxheight, ls=None, silent=False):
        """
        land use inside an area around the trajectory



        add the shape to ``self.shapes['shp_below{:.1f}km']``

        and the statistics to ``self.stat_ls['occ_shp_below{:.1f}km']``

        sizes of circle (in 24h steps):

        .. code-block:: python

           [ 0.05      ,  0.15      ,  0.4       ,  0.77666667,  1.25666667,
            1.81666667,  2.43333333,  3.08333333,  3.74333333,  4.39      ,  5.]

        Args:
            maxheight: maximum height in meters or md for MIXDEPTH
            ls (:class:`trace_source.land_sfc.land_sfc`, optional): pre loaded land surface information
                (separate loading consumes lot of time)
            silent (bool, optional): verbose output or not

        .. deprecated:: 0.1
            use the ensemble trajectories instead
        """

        import shapely.geometry as sgeom
        from shapely.ops import cascaded_union

        if maxheight == 'md':
            maxheight = self.data[1]['MIXDEPTH']
            maxheightstr = 'md'
        else:
            maxheightstr = '{:.1f}'.format(maxheight/1000.)

        lat = self.data[1]['latitude'][self.data[1]['height'] < maxheight]
        lon = self.data[1]['longitude'][self.data[1]['height'] < maxheight]
        age = self.data[1]['age'][self.data[1]['height'] < maxheight]

        if lat.shape[0] > 0:
            r = np.abs(age)**2*(5./240**2)
            p = np.poly1d([ -2.81314300e-07, 1.50462963e-04, 7.17592593e-04, 5.00000000e-02])
            r = p(abs(age))
            r[r<0.05] = 0.05

            #track = sgeom.LineString(zip(lon, lat))
            #track = track.buffer(0.1)

            shape = []
            for coord in zip(lat, lon, r):
                shape.append(sgeom.Point(coord[1], coord[0]).buffer(coord[2]))

            #shape.append(track)
            shape = cascaded_union(shape)

            # check if shape is a Polygon; if yes convert it to Multipolygon
            if shape.geom_type == 'Polygon':
                shape = sgeom.MultiPolygon([shape])



            if ls is None:
                ls = trace_source.land_sfc.land_sfc()

            occ_stat = namedtuple('occ_stat', 'no_below counter')
            cat = ls.get_land_sfc_shape(shape).compressed()
            no = float(cat.shape[0]) if cat.shape[0] > 0 else -1
            c = {x: cat.tolist().count(x) / float(no) for x in list(ls.categories.keys())}

            key = 'shp_below{}km'.format(maxheightstr)
            self.shapes[key] = shape

            key = 'occ_shp_below{}km'.format(maxheightstr)
            self.stat_ls[key] = occ_stat(no_below=no, counter=c)

            return shape

        else:
            #print('! nothin below ', maxheight)
            print('! nothin below ')
            # how do we treat this?


        
    def evaluate(self, silent=False):
        """
        adds some general statistic to a single trajectory

        Keys: ``Tmin_whole``, ``Precip_whole``, ``RHmax_whole``, ``Tmin_24h``, ``Precip_24h``,
        ``RHmax24h``, ``min_height_ag``, ``traj_distance_total``, ``mean_bearing_from_endpoint``

        Args:
            silent (bool, optional): verbose output or not

        """

        self.statistics['Tmin_whole'] = np.min(self.data[1]['AIR_TEMP'][:]-273.15)
        self.statistics['Precip_whole'] = np.sum(self.data[1]['RAINFALL'][:])
        self.statistics['RHmax_whole'] = np.max(self.data[1]['RELHUMID'][:])

        index = np.where(self.data[1]['age']==-24)[0][0]
        self.statistics['Tmin_24h'] = np.min(self.data[1]['AIR_TEMP'][:index+1]-273.15)
        self.statistics['Precip_24h'] = np.sum(self.data[1]['RAINFALL'][:index+1])
        self.statistics['RHmax_24h'] = np.max(self.data[1]['RELHUMID'][:index+1])

        self.statistics['min_height_ag'] = np.min(self.data[1]['height'])

        # distance between two points
        # bearing to start point
        dist, angle = [], []
        endcoord = (self.info['lat'], self.info['lon'])
        coords = list(zip(self.data[1]['latitude'], self.data[1]['longitude']))
        for i in range(1, len(coords)):
            d = trace_source.helpers.distance(coords[i-1], coords[i])
            a = trace_source.helpers.calculate_initial_angle(endcoord, coords[i])
            dist.append(d)
            angle.append(a)

        print('total range', sum(dist))
        self.statistics['traj_distance_total'] = sum(dist)

        comp_angle = np.exp(1j * np.array(angle))
        mean_angle = np.mean(comp_angle.real) + 1j*np.mean(comp_angle.imag)
        mean_bear = (np.angle(mean_angle, deg=True) + 360) % 360
        print('mean bear ', mean_bear)
        self.statistics['mean_bearing_from_endpoint'] = mean_bear

        print('statistics') if not silent else None
        for k, v in self.statistics.items():
            print(k, v) if not silent else None


        
if __name__ == '__main__':

    #traj = trajectory()
    #traj.load_file('../prior_examples/test_trajectories/hysplit_trajectory-20170225-00-34_677-33_038-4500_0-240.tdump')

    # traj = trajectory()
    # traj.load_file('../prior_examples/trajectories_limassol/hysplit_trajectory-20170225-03-34_677-33_038-1500_0-240.tdump')
    # traj.add_land_sfc()
    # traj.add_ensemble_land_sfc()

    traj = trajectory(config='../config_macca.toml')
    #traj.load_file('../prior_examples/trajectories_barbados/hysplit_trajectory-20140223-21-13_15--59_26-10500_0-240.tdump')
    #traj.load_file('../prior_examples/trajectories_limassol/hysplit_trajectory-20170304-12-34_677-33_038-12000_0-240.tdump')
    traj.load_file('../prior_examples/trajectories_macca/hysplit_trajectory-20160531-15--54_6199-158_861-3000_0-240.tdump')
    plot_trajectories_ens(traj, '../prior_examples/plots')
