
# coding: utf-8

import datetime, math
import os
import gc
import argparse
import numpy as np
#from Scientific.IO import NetCDF
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ast
import toml

plt.rcParams['svg.fonttype'] = 'none'

# do we need viridis?
#print(os.listdir('output/'))
#print(os.listdir('output/barbados/'))

import logging
logger = logging.getLogger(__name__)

def gen_file_list(start, end, path, station, model):
    l = []
    current = start
    while current <= end:
        #l.append('output/{}_{}_hysplit-output.nc'.format(current.strftime('%Y%m%d'), station))
        #l.append('output/{1}/{0}_{1}_hysplit-output.nc'.format(current.strftime('%Y%m%d'), station))
        l.append(path + '/{0}_{1}_{2}-output.nc'.format(current.strftime('%Y%m%d'), station, model))
        current = current + datetime.timedelta(days=1)
    return l


def plot_landsfc_2d(f, parameter, dt_list, config, savepath, config_dict, model, fcst):
    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]
    height_list = f.variables["range"][:]
    no_plots = len(dt_list)

    fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(12, 6))
    #fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(9, 6))
    ls_colors = ['lightskyblue', 'darkgreen', 'khaki', 'palegreen', 'red', 'gray', 'tan']
    ls_colors = ['lightskyblue', 'seagreen', 'khaki', '#6edd6e', 'darkmagenta', 'gray', 'tan']

    for it, dt in enumerate(dt_list):

        occ_height = f.variables[parameter][it, :, :]
        occ_height = np.ma.masked_less(occ_height, 0)
        occ_left = np.cumsum(occ_height, axis=1)

        categories = ast.literal_eval(f.variables[parameter].comment)

        # fix the typo
        if categories[1] == 'forrest':
            categories[1] == 'forest'

        l1 = axes[it].barh(height_list, occ_height[:, 0].T, 
                           align='center', height=0.35, color=ls_colors[0], edgecolor='none')
        l2 = axes[it].barh(height_list, occ_height[:, 1].T, left=occ_left[:, 0].T,
                           align='center', height=0.35, color=ls_colors[1], edgecolor='none')
        l3 = axes[it].barh(height_list, occ_height[:, 2].T, left=occ_left[:, 1].T,
                           align='center', height=0.35, color=ls_colors[2], edgecolor='none')
        l4 = axes[it].barh(height_list, occ_height[:, 3].T, left=occ_left[:, 2].T,
                           align='center', height=0.35, color=ls_colors[3], edgecolor='none')
        l5 = axes[it].barh(height_list, occ_height[:, 4].T, left=occ_left[:, 3].T,
                           align='center', height=0.35, color=ls_colors[4], edgecolor='none')
        l6 = axes[it].barh(height_list, occ_height[:, 5].T, left=occ_left[:, 4].T,
                           align='center', height=0.35, color=ls_colors[5], edgecolor='none')

        l7 = axes[it].barh(height_list, occ_height[:, 6].T, left=occ_left[:, 5].T,
                          align='center', height=0.35, color=ls_colors[6], edgecolor='none')

        #axes[it].set_ylim([0, 12])
        axes[it].set_ylim([0, config['height']['plottop']/1000.])
        axes[it].tick_params(axis='y', which='major', labelsize=13, 
                             width=1.5, length=3)
        axes[it].tick_params(axis='both', which='minor', width=1, length=2)
        axes[it].tick_params(axis='both', which='both', right=True, top=True,
                             direction='in')
        axes[it].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1.0))
        axes[it].set_xlabel(dt.strftime('%H:%M'), fontsize=13)
        
        axes[it].set_xlim(right=1)
        axes[it].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        dsp_text = 'norm. residence time'
        #axes[it].xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0, xright]))
        axes[it].tick_params(axis='x', labeltop=False, labelbottom=True)

    param_str = parameter.replace('rt_normed_', '') 
    axes[0].set_ylabel("Height [km]", fontweight='semibold', fontsize=14)
    axes[-1].tick_params(axis='x', labeltop=True, labelbottom=False, labelsize=11)
    plt.suptitle("{}   {}   {}".format(dt.strftime("%Y%m%d"), config_dict['station']["name"], param_str), 
                 fontweight='semibold', fontsize=15)
    if 'moving' in config_dict['station'].keys() and config_dict['station']['moving']:
        pass
    else:
        axes[-1].annotate("Endpoint: {:.1f} {:.1f} ".format(config_dict['station']["lat"], 
                                                            config_dict['station']["lon"]), 
                            xy=(.91, 0.96), xycoords='figure fraction',
                            horizontalalignment='center', verticalalignment='bottom',
                            fontsize=12)
    axes[-1].annotate("{} ".format(model), 
                      xy=(.91, 0.925), xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12)
    if fcst:
        axes[-1].annotate("forecast data ",
                          xy=(.7, 0.0), xycoords='figure fraction',
                          horizontalalignment='center', verticalalignment='bottom',
                          color='red', fontsize=14)
    axes[0].annotate('Time UTC', xy=(.5, .0),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=14, fontweight='semibold')
    axes[-1].annotate(dsp_text, xy=(.90, 0.86),
                      xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12, fontweight='semibold')
    fig.legend((l1, l2, l3, l4, l5, l6, l7), list(categories.values()),
               #bbox_to_anchor=(0.85, 0.952),
               bbox_to_anchor=(0.03, 0.947), loc='upper left',
               ncol=4, fontsize=11, framealpha=0.8)

    #fig.set_tight_layout({'rect': [0, 0, 1, 1], 'pad': 0.1, 'h_pad': 1.5})
    #plt.tight_layout(w_pad=0.0002)
    plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    plt.tight_layout(rect=[0, 0.02, 1, 0.88])
    fig.subplots_adjust(wspace=0, top=0.80)

    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.png".format(short_name, param_str.replace('_', '-'))
    fig.savefig(savename, dpi=400)
    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.svg".format(short_name, param_str.replace('_', '-'))
    fig.savefig(savename)
    plt.close()


def plot_regions_2d(f, parameter, dt_list, config, savepath, config_dict, model, fcst):
    with open('geonames_config.toml') as config_file:
        geo_config = toml.loads(config_file.read())

    geo_names = {int(k): v for k, v in geo_config[config['geonames']]['geo_names'].items()}
    #assert geo_names == ast.literal_eval(f.variables[parameter].comment), "geonames not matching {}".format(str(geo_names))

    NUM_COLORS = 7
    cm = plt.get_cmap('Set2')
    cNorm  = matplotlib.colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    colors = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]
    colors = [(0.65098039215686276, 0.84705882352941175, 0.32941176470588235, 1.0), 
              (1.0, 0.85098039215686272, 0.18431372549019609, 1.0), 
              (0.89803921568627454, 0.7686274509803922, 0.58039215686274515, 1.0),
              (0.40000000000000002, 0.76078431372549016, 0.6470588235294118, 1.0), 
              (0.9882352941176471, 0.55294117647058827, 0.3843137254901961, 1.0), 
              (0.55294117647058827, 0.62745098039215685, 0.79607843137254897, 1.0),  
              (0.70196078431372544, 0.70196078431372544, 0.70196078431372544, 1.0)]

    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]    
    height_list = f.variables["range"][:]
    no_plots = len(dt_list)

    fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(12, 6))
    #fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(9, 6))

    for it, dt in enumerate(dt_list):
        occ_height = f.variables[parameter][it, :, :]
        occ_height = np.ma.masked_less(occ_height, 0)
        occ_left = np.cumsum(occ_height, axis=1)

        l = []
        for i in range(len(geo_names)):
            if i == 0:
                l.append(axes[it].barh(height_list, occ_height[:, 0].T, left=0, 
                    align='center', height=0.3, color=colors[0], edgecolor='none'))
            else:
                l.append(axes[it].barh(height_list, occ_height[:, i].T, left=occ_left[:, i-1].T,
                    align='center', height=0.3, color=colors[i], edgecolor='none'))

        axes[it].set_ylim([0, 12])
        axes[it].set_ylim([0, config['height']['plottop']/1000.])
        axes[it].tick_params(axis='y', which='major', labelsize=14, 
                             width=1.5, length=3)
        axes[it].tick_params(axis='both', which='minor', width=1, length=2)
        axes[it].tick_params(axis='both', which='both', right=True, top=True,
                             direction='in')
        axes[it].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1.0))
        axes[it].set_xlabel(dt.strftime('%H:%M'), fontsize=14)
        
        axes[it].set_xlim(right=1)
        axes[it].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        dsp_text = 'norm. residence time'
        #axes[it].xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0, xright]))
        axes[it].tick_params(axis='x', labeltop=False, labelbottom=False)


    param_str = parameter.replace('rt_normed_', '') 
    axes[0].set_ylabel("Height [km]", fontweight='semibold', fontsize=14)
    axes[-1].tick_params(axis='x', labeltop=True, labelbottom=False, labelsize=11)
    plt.suptitle("{}   {}   {}".format(dt.strftime("%Y%m%d"), config_dict['station']["name"], param_str), 
                 fontweight='semibold', fontsize=15)
    if 'moving' in config_dict['station'].keys() and config_dict['station']['moving']:
        pass
    else:
        axes[-1].annotate("Endpoint: {:.1f} {:.1f} ".format(config_dict['station']["lat"], 
                                                            config_dict['station']["lon"]), 
                        xy=(.91, 0.96), xycoords='figure fraction',
                        horizontalalignment='center', verticalalignment='bottom',
                        fontsize=12)
    axes[-1].annotate("{} ".format(model), 
                      xy=(.91, 0.925), xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12)
    if fcst:
        axes[-1].annotate("forecast data ", 
                          xy=(.7, 0.0), xycoords='figure fraction',
                          horizontalalignment='center', verticalalignment='bottom',
                          color='red', fontsize=14)
    axes[0].annotate('Time UTC', xy=(.5, .0),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=14, fontweight='semibold')
    axes[-1].annotate(dsp_text, xy=(.90, 0.86),
                      xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12, fontweight='semibold')

    fig.legend(l, list(geo_names.values()),
               loc='upper left',
               #bbox_to_anchor=(0.01, 0.952),
               bbox_to_anchor=(0.01, 0.947),
               ncol=5, fontsize=12, framealpha=0.8)
    #fig.set_tight_layout({'rect': [0, 0, 1, 1], 'pad': 0.1, 'h_pad': 1.5})
    #plt.tight_layout(w_pad=0.0002)

    #plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    #plt.tight_layout(rect=[0, 0.02, 1, 0.90])
    plt.tight_layout(rect=[0, 0.02, 1, 0.88])
    fig.subplots_adjust(wspace=0, top=0.80)

    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.png".format(short_name, param_str.replace('_', '-'))
    fig.savefig(savename, dpi=400)
    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.svg".format(short_name, param_str.replace('_', '-'))
    fig.savefig(savename)
    plt.close()



def plot_lat_2d(f, parameter, dt_list, config, savepath, config_dict, model, fcst):


    #geo_names = {int(k): v for k, v in geo_config[config['geonames']]['geo_names'].items()}
    #assert geo_names == ast.literal_eval(f.variables[parameter].comment), "geonames not matching {}".format(str(geo_names))
    #lat_names = {0: 'N -60', 1: 'S -60'}

    lat_names = ast.literal_eval(f.variables[parameter].comment)

    NUM_COLORS = f.variables[parameter].shape[-1]
    #cm = plt.get_cmap('Set2')
    #cNorm  = matplotlib.colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    #scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
    #colors = [scalarMap.to_rgba(i) for i in range(NUM_COLORS)]
    #colors = [matplotlib.cm.get_cmap('Accent')(i) for i in np.linspace(0,1,NUM_COLORS)]
    colors = [
            (0.2196078431372549, 0.4235294117647059, 0.6901960784313725, 1.0), #blue
            (0.7450980392156863, 0.6823529411764706, 0.8313725490196079, 1.0), #violet
            (0.9, 0.9, 0.2, 1.0),  # yellow
            (0.7490196078431373, 0.3568627450980392, 0.09019607843137253, 1.0), #brown
            (0.4980392156862745, 0.788235294117647, 0.4980392156862745, 1.0), #green
            (0.5, 0.5, 0.5, 1.0) # grey
            ]

    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]    
    height_list = f.variables["range"][:]
    no_plots = len(dt_list)

    fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(12, 6))
    #fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(9, 6))

    for it, dt in enumerate(dt_list):
        occ_height = f.variables[parameter][it, :, :]
        occ_height = np.ma.masked_less(occ_height, 0)
        occ_left = np.cumsum(occ_height, axis=1)

        l = []
        for i in range(NUM_COLORS):
            if i == 0:
                l.append(axes[it].barh(height_list, occ_height[:, 0].T, left=0, 
                    align='center', height=0.3, color=colors[0], edgecolor='none'))
            else:
                l.append(axes[it].barh(height_list, occ_height[:, i].T, left=occ_left[:, i-1].T,
                    align='center', height=0.3, color=colors[i], edgecolor='none'))

        axes[it].set_ylim([0, 12])
        axes[it].set_ylim([0, config['height']['plottop']/1000.])
        axes[it].tick_params(axis='y', which='major', labelsize=14, 
                             width=1.5, length=3)
        axes[it].tick_params(axis='both', which='minor', width=1, length=2)
        axes[it].tick_params(axis='both', which='both', right=True, top=True,
                             direction='in')
        axes[it].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1.0))
        axes[it].set_xlabel(dt.strftime('%H:%M'), fontsize=14)
        
        axes[it].set_xlim(right=1.05)
        axes[it].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        dsp_text = 'norm. residence time'
        #axes[it].xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0, xright]))
        axes[it].tick_params(axis='x', labeltop=False, labelbottom=False)

    param_str = parameter.replace('rt_normed_', '') 
    axes[0].set_ylabel("Height [km]", fontweight='semibold', fontsize=14)
    axes[-1].tick_params(axis='x', labeltop=True, labelbottom=False, labelsize=11)
    plt.suptitle("{}   {}   {}".format(dt.strftime("%Y%m%d"), config_dict['station']["name"], param_str), 
                 fontweight='semibold', fontsize=15)
    if 'moving' in config_dict['station'].keys() and config_dict['station']['moving']:
        pass
    else:
        axes[-1].annotate("Endpoint: {:.1f} {:.1f} ".format(config_dict['station']["lat"], 
                                                            config_dict['station']["lon"]), 
                        xy=(.91, 0.96), xycoords='figure fraction',
                        horizontalalignment='center', verticalalignment='bottom',
                        fontsize=12)
    axes[-1].annotate("{} ".format(model), 
                      xy=(.91, 0.925), xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12)
    if fcst:
        axes[-1].annotate("forecast data ",
                          xy=(.7, 0.0), xycoords='figure fraction',
                          horizontalalignment='center', verticalalignment='bottom',
                          color='red', fontsize=14)
    axes[0].annotate('Time UTC', xy=(.5, .0),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=14, fontweight='semibold')
    axes[-1].annotate(dsp_text, xy=(.90, 0.86),
                      xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12, fontweight='semibold')

    fig.legend(l, list(lat_names.values()),
               loc='upper left',
               #bbox_to_anchor=(0.01, 0.952),
               bbox_to_anchor=(0.01, 0.947),
               ncol=6, fontsize=12, framealpha=0.8)
    #fig.set_tight_layout({'rect': [0, 0, 1, 1], 'pad': 0.1, 'h_pad': 1.5})
    #plt.tight_layout(w_pad=0.0002)

    #plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    #plt.tight_layout(rect=[0, 0.02, 1, 0.90])
    plt.tight_layout(rect=[0, 0.02, 1, 0.88])
    fig.subplots_adjust(wspace=0, top=0.80)

    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.png".format(short_name, param_str.replace('_', '-'))
    fig.savefig(savename, dpi=400)
    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.svg".format(short_name, param_str.replace('_', '-'))
    fig.savefig(savename)
    plt.close()



def plot_seaice_2d(f, parameter, dt_list, config, savepath, config_dict, model, fcst):

    lat_names = ast.literal_eval(f.variables[parameter].comment)

    NUM_COLORS = f.variables[parameter].shape[-1]
    # "{0: \'outofarea\', 1: \'open ocean\', 2: \'1..30\', 3: \'30..80\', 4: \'>80\', 5: \'mask\'}" 
    colors = [
            (0.41, 0.41, 0.41, 1.0), # grey
            (0.0, 0.0, 0.6, 1.0), # darkblue
            (0.2, 0.4, 0.8, 1.0), # blue
            (0.0, 0.6, 1.0, 1.0), # lightblue
            (0.8, 1.0, 1.0, 1.0), # cyan
            (0.8, 0.8, 0.8, 1.0), # grey
            ]

    time_list = f.variables["timestamp"][:]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]    
    height_list = f.variables["range"][:]
    no_plots = len(dt_list)

    fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(12, 6))
    #fig, axes = plt.subplots(1, no_plots, sharex=True, sharey=True, figsize=(9, 6))

    for it, dt in enumerate(dt_list):
        occ_height = f.variables[parameter][it, :, :]
        occ_height = np.ma.masked_less(occ_height, 0)
        occ_left = np.cumsum(occ_height, axis=1)

        l = []
        for i in range(NUM_COLORS):
            if i == 0:
                l.append(axes[it].barh(height_list, occ_height[:, 0].T, left=0, 
                    align='center', height=0.3, color=colors[0], edgecolor='none'))
            else:
                l.append(axes[it].barh(height_list, occ_height[:, i].T, left=occ_left[:, i-1].T,
                    align='center', height=0.3, color=colors[i], edgecolor='none'))

        axes[it].set_ylim([0, 12])
        axes[it].set_ylim([0, config['height']['plottop']/1000.])
        axes[it].tick_params(axis='y', which='major', labelsize=14, 
                             width=1.5, length=3)
        axes[it].tick_params(axis='both', which='minor', width=1, length=2)
        axes[it].tick_params(axis='both', which='both', right=True, top=True,
                             direction='in')
        axes[it].yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(1.0))
        axes[it].set_xlabel(dt.strftime('%H:%M'), fontsize=14)
        
        axes[it].set_xlim(right=1.05)
        axes[it].xaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
        dsp_text = 'norm. residence time'
        #axes[it].xaxis.set_minor_locator(matplotlib.ticker.FixedLocator([0, xright]))
        axes[it].tick_params(axis='x', labeltop=False, labelbottom=False)

    param_str = parameter.replace('rt_normed_', '') 
    axes[0].set_ylabel("Height [km]", fontweight='semibold', fontsize=14)
    axes[-1].tick_params(axis='x', labeltop=True, labelbottom=False, labelsize=11)
    plt.suptitle("{}   {}   {}".format(dt.strftime("%Y%m%d"), config_dict['station']["name"], param_str), 
                 fontweight='semibold', fontsize=15)
    if 'moving' in config_dict['station'].keys() and config_dict['station']['moving']:
        pass
    else:
        axes[-1].annotate("Endpoint: {:.1f} {:.1f} ".format(config_dict['station']["lat"], 
                                                            config_dict['station']["lon"]), 
                        xy=(.91, 0.96), xycoords='figure fraction',
                        horizontalalignment='center', verticalalignment='bottom',
                        fontsize=12)
    axes[-1].annotate("{} ".format(model), 
                      xy=(.91, 0.925), xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12)
    if fcst:
        axes[-1].annotate("forecast data ",
                          xy=(.7, 0.0), xycoords='figure fraction',
                          horizontalalignment='center', verticalalignment='bottom',
                          color='red', fontsize=14)
    axes[0].annotate('Time UTC', xy=(.5, .0),
                     xycoords='figure fraction',
                     horizontalalignment='center', verticalalignment='bottom',
                     fontsize=14, fontweight='semibold')
    axes[-1].annotate(dsp_text, xy=(.90, 0.86),
                      xycoords='figure fraction',
                      horizontalalignment='center', verticalalignment='bottom',
                      fontsize=12, fontweight='semibold')

    fig.legend(l, list(lat_names.values()),
               loc='upper left',
               #bbox_to_anchor=(0.01, 0.952),
               bbox_to_anchor=(0.01, 0.947),
               ncol=6, fontsize=12, framealpha=0.8)
    #fig.set_tight_layout({'rect': [0, 0, 1, 1], 'pad': 0.1, 'h_pad': 1.5})
    #plt.tight_layout(w_pad=0.0002)

    #plt.tight_layout(rect=[0, 0.02, 1, 0.93])
    #plt.tight_layout(rect=[0, 0.02, 1, 0.90])
    plt.tight_layout(rect=[0, 0.02, 1, 0.88])
    fig.subplots_adjust(wspace=0, top=0.80)

    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.png".format(short_name, param_str.replace('_', '-'))
    print(savename)
    fig.savefig(savename, dpi=400)
    savename = savepath + "/" + dt.strftime("%Y%m%d") + "_{}_{}.svg".format(short_name, param_str.replace('_', '-'))
    print(savename)
    fig.savefig(savename)
    plt.close()

def plot_filename(filename, config_dict, model, fcst, config_file='config.toml'):
    f = netCDF4.Dataset(filename, 'r')
    time_list = f.variables["timestamp"][:]
    height_list = f.variables["range"][:]

    with open(config_file) as config_file:
        config = toml.loads(config_file.read())

    #print(time_list)
    time_delta = (time_list[1:] - time_list[:-1])[0]
    dt_list = [datetime.datetime.fromtimestamp(time) for time in time_list]
    #print(dt_list)

    #savepath = 'plots/{}/'.format(dt_list[0].strftime('%Y%m%d'))
    if model != 'hysplit':
        savepath = '{}/{}_{}/'.format(
            config['plot_dir'], dt_list[0].strftime('%Y%m%d'), model)
    else:
        savepath = '{}/{}/'.format(config['plot_dir'], dt_list[0].strftime('%Y%m%d'))
    print('savepath ', savepath)
    if not os.path.isdir(savepath):
        os.makedirs(savepath)

    ## coordinate arrays for pcmesh
    #print(height_list)
    height_pcmesh = (height_list - np.mean((height_list[1:] - height_list[:-1])/2.)).tolist()
    height_pcmesh.append(height_pcmesh[-1] + np.mean((height_list[1:] - height_list[:-1])))
    #print(height_pcmesh)
    dt_pcmesh = [dt-datetime.timedelta(seconds=int(time_delta)/2) for dt in dt_list]
    dt_pcmesh.append(dt_pcmesh[-1]+datetime.timedelta(seconds=int(time_delta)))
    #print(dt_pcmesh)

    print(f)
    #print(f.variables.keys())


    if 'min_height_ag' in f.variables.keys():
        min_height_ag = f.variables['min_height_ag'][:]/1000.

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(min_height_ag),
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Trajectory min Height [km]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M") + "_{}_min_height.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'Tmin_whole' in f.variables.keys():
        var = f.variables['Tmin_whole'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=-55, vmax=5,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Trajectory min T [C]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M")+ "_{}_min_T_whole.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'Tmin_24h' in f.variables.keys():
        var = f.variables['Tmin_24h'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=-55, vmax=5,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Trajectory min T 24h [C]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M") + "_{}_mint_T_24h.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'mean_bearing_from_endpoint' in f.variables.keys():
        var = f.variables['mean_bearing_from_endpoint'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=0, vmax=360,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Bearing of Trajectory (mean)", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M")+ "_{}_mean_bearing.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()

    if 'traj_distance_total' in f.variables.keys():
        var = f.variables['traj_distance_total'][:]

        fig, ax = plt.subplots(1, figsize=(10, 5.7))
        pcmesh = ax.pcolormesh(matplotlib.dates.date2num(dt_pcmesh[:]),
                            height_pcmesh[:],
                            np.transpose(var),
                            vmin=2500, vmax=18000,
                            cmap='gist_rainbow')
        cbar = fig.colorbar(pcmesh)
        ax.set_xlim([dt_list[0], dt_list[0]+datetime.timedelta(hours=24)])
        ax.set_ylim([0, height_list[-1]])
        ax.set_xlabel("Time UTC", fontweight='semibold', fontsize=15)
        ax.set_ylabel("Height [km]", fontweight='semibold', fontsize=15)
        cbar.ax.set_ylabel("Total length [km]", fontweight='semibold', fontsize=15)
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))

        ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=3))
        ax.xaxis.set_minor_locator(matplotlib.dates.MinuteLocator(byminute=[0,30]))
        ax.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(500))

        ax.tick_params(axis='both', which='major', labelsize=14, 
                    width=2, length=4)
        ax.tick_params(axis='both', which='minor', width=1.5, length=3)
        ax.tick_params(axis='both', which='both', right=True, top=True,
                    direction='in')
        cbar.ax.tick_params(axis='both', which='major', labelsize=14,
                            width=2, length=4)

        savename = savepath + "/" + dt_list[0].strftime("%Y%m%d_%H%M") + "_{}_length_tot.png".format(short_name)
        fig.savefig(savename, dpi=250)
        plt.close()


    for rh in config['height']['reception']:
        if rh != 'md':
            rh_string = rh + 'km'
        else:
            rh_string = rh
        
        if 'rt_normed_lat_below'+rh_string in f.variables:
            plot_lat_2d(f, 'rt_normed_lat_below'+rh_string, dt_list, config, savepath, config_dict, model, fcst)
        plot_regions_2d(f, 'rt_normed_region_below'+rh_string, dt_list, config, savepath, config_dict, model, fcst)
        plot_landsfc_2d(f, 'rt_normed_landsfc_below'+rh_string, dt_list, config, savepath, config_dict, model, fcst)
        if 'rt_normed_seaice_below'+rh_string in f.variables:
            plot_seaice_2d(f, 'rt_normed_seaice_below'+rh_string, dt_list, config, savepath, config_dict, model, fcst)


    gc.collect()

parser = argparse.ArgumentParser(description='Plot trace output')
parser.add_argument('--station', help='station name like limassol, barbados or mcmurdo', required=True)
parser.add_argument('--date', help='date in the format YYYYMMDD')
parser.add_argument('--daterange', help='date range in the format YYYYMMDD-YYYMMDD')
parser.add_argument('--model', help='model to plot, default hysplit', default='hysplit')
parser.add_argument('--fcst', default='false', help='add a flag indicating forecast data')

args = parser.parse_args()

if args.fcst == 'false':
    fcst = False
elif args.fcst == 'true':
    fcst = True
else:
    raise ValueError

if args.date is not None:
    begin = datetime.datetime.strptime(args.date, '%Y%m%d')
    end = begin
if args.daterange is not None:
    begin = datetime.datetime.strptime(args.daterange.split('-')[0], '%Y%m%d')
    end = datetime.datetime.strptime(args.daterange.split('-')[1], '%Y%m%d')
station = args.station

if args.station is not None:
    config = 'config_{}.toml'.format(args.station)

config_dict = toml.load(config)
short_name = config_dict['station']['short_name']

#filelist = gen_file_list(datetime.date(2016, 1, 22), datetime.date(2016, 1, 24), station)
filelist = gen_file_list(begin, end, config_dict['output_dir'], short_name, args.model)
for filename in filelist:
    print('=============================================================================')
    print('plotting file ', filename)
    plot_filename(filename, config_dict, args.model, fcst, config_file=config)

logger.warning('DeprecationWarning: New version of plot2d.py. For old version see plot2d_legacy.py')
