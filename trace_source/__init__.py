#! /usr/bin/env python3
# coding=utf-8
""""""
"""
Author: radenz@tropos.de

"""

import datetime
import toml
import numpy as np
import netCDF4
import os
import gc
import subprocess



def time_list(begin, end, delta):
    """generate a time list from begin to <= end with delta

    Args:
        begin (datetime): start time
        end (datetime): end time (included)

    Returns:
        list: List of generated datetime objects

    """
    time_list = []
    elem = begin
    while elem <= end:
        time_list.append(elem)
        elem += datetime.timedelta(hours=delta)
    return time_list


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

def get_git_hash():
    """
    Returns:
        git describe string
    """
    commit_id = subprocess.check_output(['git', 'describe', '--always'])
    branch = subprocess.check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
    return commit_id.rstrip(), branch.rstrip()


class assemble_pattern():
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

        self.no_part = []
        self.time_res = []

    # the assemble method is model data specific
    #def assemble(self, dt_range=None):


    def dump2netcdf(self, model_str='hysplit'):
        """
        dump the assembled data to a netcdf file repeatatly call :meth:`save_item`
        configuration (directories, names, etc) is given by the config file
        """

        timestamps = np.array([(dt - datetime.datetime(1970,1,1)).total_seconds() for dt in self.dt_list])
        hours_cn = np.array([dt.hour + dt.minute / 60. + dt.second / 3600. for dt in self.dt_list])

        if not os.path.isdir(self.config['output_dir']):
            os.makedirs(self.config['output_dir'])

        ncfile = self.config['output_dir'] +\
                '{}_{}_{}-output.nc'.format(
                    self.dt_list[0].strftime('%Y%m%d'), self.config['station']['short_name'], model_str)
        # ncfile = "/home/devel/" +\
        #          '{}_hysplit_output.nc'.format(self.config['time']['begin_dt'].strftime('%Y%m%d'))

        #dataset = netCDF4.Dataset(ncfile, 'w', format='NETCDF4')
        dataset = netCDF4.Dataset(ncfile, 'w', format='NETCDF3_CLASSIC')

        dim_time = dataset.createDimension('time', len(self.dt_list))
        dim_height = dataset.createDimension('height', len(self.height_list))
        dim_age = dataset.createDimension('time_age', abs(self.config['time']['tr_duration'])+1)
        dim_cat = dataset.createDimension('categories', 7)
        dim_regions = dataset.createDimension('regions', len(list(self.geo_names.keys())))
        dim_lats = dataset.createDimension('lat_thres', len(list(self.lat_names.keys())))

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

        save_item(dataset, {'var_name': 'no_part', 'dimension': ('time', ),
                            'arr': np.array(self.no_part),
                            'long_name': "number particles/trajectories",
                            'units': "no"})
        save_item(dataset, {'var_name': 'time_res', 'dimension': ('time', ),
                            'arr': np.array(self.time_res),
                            'long_name': "backward time resolution",
                            'units': "no"})

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
                                 'long_name': "land surface " + key.lower(),
                                 'comment': str(self.ls_categories)} for key in ls_data_keys}
        for k in ls_data_keys:
            print(k, self.statls_dict.get(k).shape)
            dataset = save_item(dataset, {'var_name': modified_params[k]['var_name'],
                                          'dimension': ('time', 'height', 'categories'),
                                          'arr': self.statls_dict.get(k),
                                          'long_name': modified_params[k]['long_name'],
                                          'comment': modified_params[k]['comment']})

        for k in [ky for ky in ls_data_keys if 'ens' in ky]:
            rel = self.statls_dict.get(k)
            no_below = self.stat2d_dict.get(k + "_no_below")
            no_below = np.repeat(no_below[:,:,np.newaxis], rel.shape[-1], axis=2)
            no_below[no_below < 0] = np.nan
            norm = np.array(self.no_part)*10*(24./np.array(self.time_res))
            norm = np.repeat(norm[:,np.newaxis], rel.shape[1], axis=1)
            norm = np.repeat(norm[:,:,np.newaxis], rel.shape[2], axis=2)

            normed_time = rel*no_below/norm
            normed_time[~np.isfinite(normed_time)] = -1
            str_below = modified_params[k]['var_name'].replace("occ_ens_", "")
            var_name = "rt_normed_landsfc_" + str_below
            long_name = "normed residence time land surface " + str_below
            dataset = save_item(dataset, {'var_name': var_name,
                                          'dimension': ('time', 'height', 'categories'),
                                          'arr': normed_time,
                                          'long_name': long_name,
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
            
        for k in gn_data_keys:
            rel = self.statgn_dict.get(k)
            no_below = self.stat2d_dict.get(k + "_no_below")
            no_below = np.repeat(no_below[:,:,np.newaxis], rel.shape[-1], axis=2)
            no_below[no_below < 0] = np.nan
            norm = np.array(self.no_part)*10*(24./np.array(self.time_res))
            norm = np.repeat(norm[:,np.newaxis], rel.shape[1], axis=1)
            norm = np.repeat(norm[:,:,np.newaxis], rel.shape[2], axis=2)

            normed_time = rel*no_below/norm
            normed_time[~np.isfinite(normed_time)] = -1
            str_below = modified_params[k]['var_name'].replace("region_ens_", "")
            var_name = "rt_normed_region_" + str_below
            long_name = "normed residence time named region " + str_below
            dataset = save_item(dataset, {'var_name': var_name,
                                          'dimension': ('time', 'height', 'regions'),
                                          'arr': normed_time,
                                          'long_name': long_name,
                                          'comment': modified_params[k]['comment']})


        # TODO make statlat optional statlat_dict
        lat_data_keys = list(self.statlat_dict.keys())
        print('lat_data_keys ', lat_data_keys)
        modified_params = {key: {'var_name': key,
                                 'long_name': "lat_thres " + key.lower(),
                                 'comment': str(self.lat_names)} for key in lat_data_keys}
        for k in lat_data_keys:
            print("self.statlat_dict.keys()", self.statlat_dict.keys())
            print(k, self.statlat_dict.get(k).shape)
            dataset = save_item(dataset, {'var_name': modified_params[k]['var_name'],
                                          'dimension': ('time', 'height', 'lat_thres'),
                                          'arr': self.statlat_dict.get(k),
                                          'long_name': modified_params[k]['long_name'],
                                          'comment': modified_params[k]['comment']})

        for k in lat_data_keys:
            rel = self.statlat_dict.get(k)
            no_below = self.stat2d_dict.get(k + "_no_below")
            no_below = np.repeat(no_below[:,:,np.newaxis], rel.shape[-1], axis=2)
            no_below[no_below < 0] = np.nan
            norm = np.array(self.no_part)*10*(24./np.array(self.time_res))
            norm = np.repeat(norm[:,np.newaxis], rel.shape[1], axis=1)
            norm = np.repeat(norm[:,:,np.newaxis], rel.shape[2], axis=2)

            normed_time = rel*no_below/norm
            normed_time[~np.isfinite(normed_time)] = -1
            str_below = modified_params[k]['var_name'].replace("lat_ens_", "")
            var_name = "rt_normed_lat_" + str_below
            long_name = "normed residence time latitude " + str_below
            dataset = save_item(dataset, {'var_name': var_name,
                                          'dimension': ('time', 'height', 'lat_thres'),
                                          'arr': normed_time,
                                          'long_name': long_name,
                                          'comment': modified_params[k]['comment']})


        # save_item(dataset, {'var_name': 'width', 'dimension': ('time', 'height'),
        #                     'arr': .corr_width_reg[:].filled(), 'long_name': "Spectral width",
        #                     'comment': "Wind profiler spectral width (standard deviation) corrected by cloud radar (only Bragg contribution)",
        #                     'units': "m s-1", 'units_html': "m s<sup>-1</sup>",
        #                     'missing_value': -99., 'plot_range': (0.01, 4.),
        #                     'plot_scale': "logarithmic"})


        with open('output_meta.toml') as output_meta:
            meta_info = toml.loads(output_meta.read())

        dataset.description = meta_info['description'][model_str]
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
        dataset.git_commit, dataset.git_branch = get_git_hash()
        dataset.close()
        gc.collect()



#
# for some wired reson the submodule imports have to be done afterwards
#
import trace_source.hysplit
import trace_source.flexpart
import trace_source.land_sfc
import trace_source.helpers
import trace_source.simulations
