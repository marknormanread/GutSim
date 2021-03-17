"""
Created on 01/05/2014

Mark N. Read, 2018
"""
import gutsim.bacteria
import numpy as np
import pandas
import math
import matplotlib
matplotlib.use('Agg')  # allows matplotlib to be run on a machine without a display.
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import os
import seaborn as sns
# Used in writing data as XML files
from xml.etree.ElementTree import Element, SubElement
from xml.etree import ElementTree
from xml.dom import minidom
import time
from gutsim.stackedBarGraph import StackedBarGrapher

# For pretty plotting
sns.set(style='darkgrid', context='talk')


guild_colours = {
    'Grf': '#DB5F57',
    'Grm': '#D3DA57',
    'Gpf': '#55B453',
    'Gpm': '#57BFD2',
    'Gmf': '#5C54A5',
    'Gmm': '#DD56A8',
    'Total': '#000000'
}


def setBoxColors(bp):
    """
    function for setting the colors of the box plots pairs. This code was taken and modified from the internet on
    2014/04/30, from:
    http://stackoverflow.com/questions/16592222/matplotlib-group-boxplots.
    The box plot must be constructed in pairs, and `bp` is the return of the `boxplot` function. E.g.,
    # First boxplot pair
    bp = boxplot(A, positions = [1, 2], widths = 0.6)
    setBoxColors(bp)

    # Second boxplot pair
    bp = boxplot(B, positions = [4, 5], widths = 0.6)
    setBoxColors(bp)

    # Third boxplot pair
    bp = boxplot(C, positions = [7, 8], widths = 0.6)
    setBoxColors(bp)
    """
    try:
        if bp['boxes']:
            plt.setp(bp['boxes'][0], color='blue')
            plt.setp(bp['caps'][0], color='blue')
            plt.setp(bp['caps'][1], color='blue')
            plt.setp(bp['whiskers'][0], color='blue')
            plt.setp(bp['whiskers'][1], color='blue')
            # plt.setp(bp['fliers'][0], color='blue')  # Raises error. Not called 'fliers' anymore?
            # plt.setp(bp['fliers'][1], color='blue')
            plt.setp(bp['medians'][0], color='blue')

            plt.setp(bp['boxes'][1], color='red')
            plt.setp(bp['caps'][2], color='red')
            plt.setp(bp['caps'][3], color='red')
            plt.setp(bp['whiskers'][2], color='red')
            plt.setp(bp['whiskers'][3], color='red')
            # plt.setp(bp['fliers'][2], color='red')
            # plt.setp(bp['fliers'][3], color='red')
            plt.setp(bp['medians'][1], color='red')
    except IndexError as err:
        print("Error in changing box plot colours. No box present to change colour on?")
        print(err)

class Vessel(object):
    """
    Superclass of all bacteria-containing bio-reactors.
    """
    def __init__(self, ID=0):
        self._bacteria = []  # Stores all the currently active bacteria in the vessel.
        # Stores the quantities of nutrient available in the vessel.
        self._nutrientStores = []
        self._ID = ID

    def get_id(self):
        return self._ID


class SpatioTemporalHeatmapLogger(object):
    def __init__(self):        
        # Data not all available upfront. 
        # Most efficient way of building up the data structure needed for making the heatmap is to store would-be
        # columns as named lists in a dictionary, appending to each list as data become available. 
        # The DataFrame can then be created and input at the 'write' stage. 
        # This will be stored in a long-form format, which is then later 'pivoted' into a 2D table for heatmap ploting. 
        self.data = dict()
        # Data series to be recorded are initialised here. 
        # Guilds relative abundances (proportions). 
        self.data['Grf_prop'] = list()
        self.data['Grm_prop'] = list()
        self.data['Gpf_prop'] = list()
        self.data['Gpm_prop'] = list()
        self.data['Gmf_prop'] = list()
        self.data['Gmm_prop'] = list()
        # Guild carbon vs nitrogen limitation. 
        self.data['Grf_nLim'] = list()
        self.data['Grm_nLim'] = list()
        self.data['Gpf_nLim'] = list()
        self.data['Gpm_nLim'] = list()
        self.data['Gmf_nLim'] = list()
        self.data['Gmm_nLim'] = list()
        # Guild doubling times
        self.data['Grf_doublingTime'] = list()
        self.data['Grm_doublingTime'] = list()
        self.data['Gpf_doublingTime'] = list()
        self.data['Gpm_doublingTime'] = list()
        self.data['Gmf_doublingTime'] = list()
        self.data['Gmm_doublingTime'] = list()
        # Nutrients.
        self.data['Cr'] = list()
        self.data['Cp'] = list()
        self.data['Cm'] = list()
        self.data['Ci'] = list()
        self.data['Nf'] = list()
        self.data['Nm'] = list()
        self.data['non_fermentable'] = list()
        self.data['other_organic_matter'] = list()
        # Record the time stamp as it should appear in plots.
        self.data['timestamp'] = list()
        self.data['ns_location'] = list()  # Index the location of the nutrient store at the current timestamp

    def log(self, timestamp, nutrientStores, bacteria): 
        # Each bacterium is assigned a nutrient store from which to extract nutrients. 
        # There are likely multiple bacteria feeding from each store.
        b_per_ns = float(len(bacteria)) / len(nutrientStores)       
        for ns_index, ns in enumerate(nutrientStores):
            self.data['timestamp'].append(timestamp)
            self.data['ns_location'].append(ns_index)
            # Log nutritional information. 
            self.data['Cr'].append(ns.Cw)
            self.data['Cp'].append(ns.Cd)
            self.data['Cm'].append(ns.Cm)
            self.data['Ci'].append(ns.Ci)
            self.data['Nf'].append(ns.Nf)
            self.data['Nm'].append(ns.Nm)
            self.data['non_fermentable'].append(ns.non_fermentable)
            self.data['other_organic_matter'].append(ns.other_organic_matter)
            # Log bacterial information. 
            # Slice to extract bacteria associated with the present nutrient store.
            ns_bacteria_start_index = int(ns_index * b_per_ns)
            ns_bacteria_end_index = int((ns_index + 1) * b_per_ns) - 1
            ns_bacteria = bacteria[ns_bacteria_start_index : ns_bacteria_end_index]
            Grf = []
            Grm = []
            Gpf = []
            Gpm = []
            Gmf = []
            Gmm = []            
            for b in ns_bacteria:  # Sort bacteria by guild
                if   isinstance(b, gutsim.bacteria.Guild_rf_f) : Grf.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_rm_m) : Grm.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_pf_f) : Gpf.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_pm_m) : Gpm.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_m_m)  : Gmm.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_mf_f) : Gmf.append(b)
            # Proportional makeup of each nutrient store's bacteria. 
            if len(ns_bacteria) > 0:  # Shouldn't happen, but safety just in case. 
                Grf_prop = 100. * len(Grf) / len(ns_bacteria) 
                Grm_prop = 100. * len(Grm) / len(ns_bacteria)
                Gpf_prop = 100. * len(Gpf) / len(ns_bacteria)
                Gpm_prop = 100. * len(Gpm) / len(ns_bacteria)
                Gmf_prop = 100. * len(Gmf) / len(ns_bacteria)
                Gmm_prop = 100. * len(Gmm) / len(ns_bacteria)
            else:
                Grf_prop = float('NaN')
                Grm_prop = float('NaN')
                Gpf_prop = float('NaN')
                Gpm_prop = float('NaN')
                Gmf_prop = float('NaN')
                Gmm_prop = float('NaN')
            self.data['Grf_prop'].append(Grf_prop)
            self.data['Grm_prop'].append(Grm_prop)
            self.data['Gpf_prop'].append(Gpf_prop)
            self.data['Gpm_prop'].append(Gpm_prop)
            self.data['Gmf_prop'].append(Gmf_prop)
            self.data['Gmm_prop'].append(Gmm_prop)
            # Proportion of each guild that is nitrogen limited. 
            Grf_state  = GuildState(Grf)
            Grm_state = GuildState(Grm)
            Gpf_state = GuildState(Gpf)
            Gpm_state = GuildState(Gpm)
            Gmf_state = GuildState(Gmf)
            Gmm_state = GuildState(Gmm)
            self.data['Grf_nLim'].append(Grf_state.n_limited)
            self.data['Grm_nLim'].append(Grm_state.n_limited)
            self.data['Gpf_nLim'].append(Gpf_state.n_limited)
            self.data['Gpm_nLim'].append(Gpm_state.n_limited)
            self.data['Gmf_nLim'].append(Gmf_state.n_limited)
            self.data['Gmm_nLim'].append(Gmm_state.n_limited)
            # Doubling time of each guild
            self.data['Grf_doublingTime'].append(Grf_state.median_doubling_time)
            self.data['Grm_doublingTime'].append(Grm_state.median_doubling_time)
            self.data['Gpf_doublingTime'].append(Gpf_state.median_doubling_time)
            self.data['Gpm_doublingTime'].append(Gpm_state.median_doubling_time)
            self.data['Gmf_doublingTime'].append(Gmf_state.median_doubling_time)
            self.data['Gmm_doublingTime'].append(Gmm_state.median_doubling_time)

    def write(self, filename_prefix):
        def generate_heatmap(df, graph_fp, centre=None, cmap=None, vmin=None, vmax=None, 
                             lognorm=False, cbar_kws=None,
                             zero_treatment='1e-2'):
            """ To plot a log-transformed heatmap, specify cbar_kws AND lognorm=True. """
            if cbar_kws is None:
                cbar_kws = dict()
            else:
                cbar_kws = dict(cbar_kws)  # Copy, so changed made here do not persist between method invocations.
            if 'ticks' not in cbar_kws.keys():
                cbar_kws['ticks'] = []

            df.to_csv(graph_fp + '.csv', index=False)
            plt.clf()
            norm = None
            min_val = df.min().min()
            max_val = df.max().max()
            # Some data series (like nutrient concentrations) are best visualised on a log-scale. 
            # This has proven difficult to do, as every option carries issues.
            # I avoided matplotlib pcolormesh as it's too fiddly; rather sns.heatmap is more usable. 
            # However, sns.heatmap doesn't handle log-transforming of data very well. 
            # Tried a few things, reported here to avoid anyone else repeating these (fruitless) investigations. 
            # Ultimate solution is to use log transform, but set values=0 to 1e-2 times less than the smallest non-zero
            # values. Hence, NaN values where no nutrient stores exist (at that point in time) remain greyed out; and
            # the two-log-decade difference between zero (post-transform) and other values makes it clear where the
            # zeroes were, whilst still providing an outline of nutrient stores over time. 
            # Tried transforming the dataframe before using heatmap, which does work (e.g. using zero-friendly 
            # arcsinh), but the colourbar itself then becomes uninterpretable.
            if lognorm and max_val != min_val:
                df_zero_as_nan = df.applymap(lambda x: x if x > 0. else float('nan'))
                non_zero_min = df_zero_as_nan.min().min()
                if zero_treatment == '1e-2':
                    # Applies a quotient to avoid log(0). Alternative is to set zero to Nan, and hence not plot it. 
                    # Set zero values to two logs lower than lowest non-zero value. 
                    df = df.applymap(lambda x: non_zero_min * 1e-2 if x == 0. else x)
                    
                elif zero_treatment == 'nan':
                    # This works, but can set large swathes of the heatmap to "unvisible" to the point that the graph
                    # (and outline of the nutrient store shape) are no longer legable. 
                    df = df_zero_as_nan
                elif zero_treatment == 'quotient':
                    # Uniform quotient applied to all; however this should be data-series specific (otherwise 
                    # data-occupying region of colour bar into a minority. 
                    quotient = 1e-10
                    df = df.applymap(lambda x: x + quotient) 
                    non_zero_min = min(non_zero_min, quotient)  # Update if quotient now smallest value. 

                if df.min().min() != df.max().max():  # Comparison of (possibly) adjusted values.
                    # Log transformed space doesn't handle zero or flat (e.g. if zero & quotient applied) data well. 
                    norm = LogNorm(vmin=non_zero_min, vmax=max_val)  # Transformation of colour bar into log space. 
                    
                    # (Manually) set location of ticks for the colourbar. 
                    # At time of writing, seaborn heatmap does not handle log-tranformed data well. 
                    if len(cbar_kws['ticks']) == 0:  # Calculate if not supplied.
                        cbar_kws['ticks'] = [math.pow(10, i) for i in range(math.floor(math.log10(non_zero_min)),
                                                                     1 + math.ceil(math.log10(df.max().max())))]

            # Here for bug testing, not essential code. 
            # df.to_csv(graph_fp + '_plotted.csv', index=False)
            if max_val == min_val and vmax is None and vmin is None:  # Don't override user-set values. 
                # Colour bar messes up in this case. Give range manually. 
                cbar_kws['ticks'] = [math.floor(min_val), math.ceil(min_val)]

            if len(cbar_kws['ticks']) == 0:
                if vmax is None and vmin is None:  # Not able to assign above.
                    cbar_kws['ticks'] = [vmax, vmin]
                else:
                    cbar_kws['ticks'] = [min_val, min_val]

            print('SELECTED TICKS for ' + graph_fp)
            print(cbar_kws['ticks'])

            sns.heatmap(data=df, cbar=True, xticklabels=False, yticklabels=False, cmap=cmap, 
                        center=centre, vmin=vmin, vmax=vmax, rasterized=True, cbar_kws=cbar_kws, norm=norm)
            plt.xlabel('')
            plt.ylabel('')
            plt.savefig(graph_fp + '.png', bbox_inches='tight', dpi=300)
            # plt.savefig(graph_fp + '.eps', bbox_inches='tight')
            plt.close()
            

        diverging_cmap = LinearSegmentedColormap.from_list('c_n_diverge', ['red', 'white', 'green'], N=256)
        # For subsequent processing into a 2D grid from which heatmaps can be plotted.     
        data_df = pandas.DataFrame.from_dict(self.data, orient='columns')
        # Nutrients.
        cbar_kws = {'format': '%.0e'}        
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Cr')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_Cr', cmap='Reds', cbar_kws=cbar_kws, lognorm=True)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Cp')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_Cp', cmap='Oranges', cbar_kws=cbar_kws, lognorm=True)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Cm')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_Cm', cmap='Blues', cbar_kws=cbar_kws, lognorm=True)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Ci')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_Ci', cmap='PuRd', cbar_kws=cbar_kws, lognorm=True)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Nf')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_Nf', cmap='Purples', cbar_kws=cbar_kws, lognorm=True)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Nm')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_Nm', cmap='Greens', cbar_kws=cbar_kws, lognorm=True)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='non_fermentable')
        mx = df.max().max()
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_non_fermentable', cmap='Greys',
                         cbar_kws={
                             'format': '%.0e',
                             'ticks': list(np.array([1.0, 0.66, 0.33, 0.0]) * mx)
                         }, lognorm=False, vmin=0., vmax=mx)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='other_organic_matter')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nutrient_other_organic_matter', cmap='YlGnBu', cbar_kws=cbar_kws, lognorm=True)
        
        # Relative abundance of community.
        # Relative abundance in log scale. 
        cbar_kws = {
            'format': '%3g',
            'ticks': [100, 50, 20, 10, 5, 1]
        }
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Grf_prop')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Grf', vmin=min(cbar_kws['ticks']), vmax=100,
                         lognorm=True, cbar_kws=cbar_kws, zero_treatment='nan')
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Grm_prop')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Grm', vmin=min(cbar_kws['ticks']), vmax=100,
                         lognorm=True, cbar_kws=cbar_kws, zero_treatment='nan')
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpf_prop')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gpf', vmin=min(cbar_kws['ticks']), vmax=100,
                         lognorm=True, cbar_kws=cbar_kws, zero_treatment='nan')
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpm_prop')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gpm', vmin=min(cbar_kws['ticks']), vmax=100,
                         lognorm=True, cbar_kws=cbar_kws, zero_treatment='nan')
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmf_prop')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gmf', vmin=min(cbar_kws['ticks']), vmax=100,
                         lognorm=True, cbar_kws=cbar_kws, zero_treatment='nan')
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmm_prop')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gmm', vmin=min(cbar_kws['ticks']), vmax=100,
                         lognorm=True, cbar_kws=cbar_kws, zero_treatment='nan')
        # Relative abundance in linear scale
        # df = data_df.pivot(index='timestamp', columns='ns_location', values='Grf_prop')
        # generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Grf', vmin=0, vmax=100)
        # df = data_df.pivot(index='timestamp', columns='ns_location', values='Grm_prop')
        # generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Grm', vmin=0, vmax=100)
        # df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpf_prop')
        # generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gpf', vmin=0, vmax=100)
        # df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpm_prop')
        # generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gpm', vmin=0, vmax=100)
        # df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmf_prop')
        # generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gmf', vmin=0, vmax=100)
        # df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmm_prop')
        # generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_prop_Gmm', vmin=0, vmax=100)
        
        # Nitrogen limitation of community.
        cbar_kws['ticks'] = [100, 80, 60, 40, 20, 0]
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Grf_nLim')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nLim_Grf', cbar_kws=cbar_kws, centre=50, vmin=0, vmax=100, cmap=diverging_cmap)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Grm_nLim')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nLim_Grm', cbar_kws=cbar_kws, centre=50, vmin=0, vmax=100, cmap=diverging_cmap)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpf_nLim')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nLim_Gpf', cbar_kws=cbar_kws, centre=50, vmin=0, vmax=100, cmap=diverging_cmap)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpm_nLim')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nLim_Gpm', cbar_kws=cbar_kws, centre=50, vmin=0, vmax=100, cmap=diverging_cmap)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmf_nLim')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nLim_Gmf', cbar_kws=cbar_kws, centre=50, vmin=0, vmax=100, cmap=diverging_cmap)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmm_nLim')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_nLim_Gmm', cbar_kws=cbar_kws, centre=50, vmin=0, vmax=100, cmap=diverging_cmap)
        # Doubling times of the community. 
        # Find maximum; plot all graphs on same scale. First max works by index, second across the resulting indexes.
        max_doubling_time = data_df.loc[:, ['Grf_doublingTime', 'Grm_doublingTime', 
                                            'Gpf_doublingTime', 'Gpm_doublingTime',
                                            'Gmf_doublingTime', 'Gmm_doublingTime']].max(axis=0).max(axis=0)
        cbar_kws = {
            # 'format': '%3g',
            'ticks': [80, 60, 40, 20, 0]
        }
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Grf_doublingTime')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_doublingTime_Grf', cbar_kws=cbar_kws, cmap='YlOrBr', vmin=0, vmax=max_doubling_time)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Grm_doublingTime')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_doublingTime_Grm', cbar_kws=cbar_kws, cmap='YlOrBr', vmin=0, vmax=max_doubling_time)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpf_doublingTime')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_doublingTime_Gpf', cbar_kws=cbar_kws, cmap='YlOrBr', vmin=0, vmax=max_doubling_time)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gpm_doublingTime')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_doublingTime_Gpm', cbar_kws=cbar_kws, cmap='YlOrBr', vmin=0, vmax=max_doubling_time)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmf_doublingTime')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_doublingTime_Gmf', cbar_kws=cbar_kws, cmap='YlOrBr', vmin=0, vmax=max_doubling_time)
        df = data_df.pivot(index='timestamp', columns='ns_location', values='Gmm_doublingTime')
        generate_heatmap(df, graph_fp=filename_prefix + '/spatiotemporal_doublingTime_Gmm', cbar_kws=cbar_kws, cmap='YlOrBr', vmin=0, vmax=max_doubling_time)


class GuildState:
    def __init__(self, bacter):
        """
        Supplied with a list of Bacterium objects, this init method will sort through them and calculate what
        proportion are in each possible state.
        """
        dead = resistant = stressed = exponential = 0  # Sets all these parameters to 0.
        n_lim = c_lim = 0  # Numbers of bacteria which are either carbon or nitrogen limited.
        growth_rate = 0
        self.growth_rates = []
        self.doubling_times = []
        for b in bacter:
            if b.is_dead(): dead += 1
            elif b.is_resistant(): resistant += 1
            elif b.is_stressed(): stressed += 1
            elif b.is_exponential(): exponential += 1
            # Record state of alive bacteria.
            if not b.is_dead():
                if b.carbon_limiting() : c_lim += 1
                else : n_lim += 1
                growth_rate += b.divide_rate
                self.growth_rates.append(b.divide_rate)
                if b.divide_rate != 0.0:
                    self.doubling_times.append(1.0 / b.divide_rate)

        total = dead + resistant + stressed + exponential
        total_alive = resistant + stressed + exponential
        if total != 0:
            self.prop_dead = 100 * float(dead) / total
            self.prop_resistant = 100 * float(resistant) / total
            self.prop_stressed = 100 * float(stressed) / total
            self.prop_exponential = 100 * float(exponential) / total
        else:
            self.prop_dead = 0
            self.prop_resistant = 0
            self.prop_stressed = 0
            self.prop_exponential = 0
        
        self.n_limited = 0
        self.c_limited = 0
        self.median_growth_rate = 0
        self.median_doubling_time = 0
        self.doubling_time = 0.0
        if total_alive != 0:
            self.n_limited = 100.0 * float(n_lim) / total_alive
            self.c_limited = 100.0 * float(c_lim) / total_alive
            # The doubling time is 1 over the mean probability of division in an hour. 
            # Gives the doubling time in hours.
            if (float(growth_rate) / total_alive) != 0:
                self.doubling_time = 1.0 / (float(growth_rate) / total_alive)
            self.median_growth_rate = np.median(self.growth_rates)
            self.median_doubling_time = np.median(self.doubling_times)

    def stack_dead(self):
        return self.prop_dead

    def stack_resistant(self):
        return self.prop_dead + self.prop_resistant

    def stack_stressed(self):
        return self.prop_dead + self.prop_resistant + self.prop_stressed

    def stack_exponential(self):
        return self.prop_dead + self.prop_resistant + self.prop_stressed + self.prop_exponential

    def stack_c_limited(self):
        return self.c_limited

    def stack_n_limited(self):
        return self.c_limited + self.n_limited


class Logger:
    def __init__(self, vessel, prefix='', logPeriod=10, bacteria_graphing_start_time=float('inf'), graphPath=None):
        """
        """
        self._vessel = vessel
        self._start_time = time.time()
        # These lists are timeseries data, and each element contains only a number, corresponding to the number
        # of each particular kind of bacteria were present at that particular sample point.
        # This is all bugs of each guild-type, dead or alive. 
        self._numGrf = []
        self._numGrm = []
        self._numGpf = []
        self._numGpm = []
        self._numGmm = []
        self._numGmf = []
        self._totalMicrobes = []
        self.Cw = []; self.Cd = []; self.Cm = []; self.Nf = []; self.Nm = []
        # Arrays for holding time series data on proportion of each guild that is in each possible state.
        self._Grf_state = []
        self._Grm_state = []
        self._Gpf_state = []
        self._Gpm_state = []
        self._Gmm_state = []
        self._Gmf_state = []
        self._time = []
        self._logCounter = 0
        self._logPeriod = logPeriod
        self._prefix = prefix  # Prefix used for file names. E.g. can be a directory in which to place data.
        self.bacteria_graphing_start_time = bacteria_graphing_start_time
        self._set_graph_path(graphPath)

        self.spatialNutrientNf = []
        self.spatialNutrientNm = []
        self.spatialNutrientCw = []
        self.spatialNutrientCd = []
        self.spatialNutrientCm = []
        # Will store a timeseries list of Nutrient objects, detailing quantities input every time step.
        self.nutrient_inputs = []

        self.heatmap_logger = SpatioTemporalHeatmapLogger()

    def _set_graph_path(self, path):
        """
        Sets the path to which vessel-specific time-snap-shot data will be written.
        """
        if not os.path.isdir(path):
            os.makedirs(path)
        self._graphPath = path

    def log(self, bacter, time, nutrient_input=None, bacteria_excreted=None):
        """
        Logs the state of the mouse. This can be done periodically, rather than every time step, since there are a
        great many timesteps and logging that frequently can slow the simulation substantially.
        Arguments:
          bacter - list containing Bacteria objects.
          time - the current simulation time.
          graph_global_dynamics - when true, print graphs of current status to the filesystem.
          bacteria_excreted (int) number of bacteria excreted since last logging event.
        """
        # Don't log every timestep, just some of them. Huge amount of data otherwise.
        if (self._logCounter % self._logPeriod) == 0:
            Grf = [];  Gpf = [];  Grm = [];  Gpm = [];  Gmm = [];  Gmf = []
            # Sort bacteria by guild
            for b in bacter:  # Count each of the bacteria, and catalogue by class.
                if   isinstance(b, gutsim.bacteria.Guild_rf_f) : Grf.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_rm_m) : Grm.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_pf_f) : Gpf.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_pm_m) : Gpm.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_m_m) : Gmm.append(b)
                elif isinstance(b, gutsim.bacteria.Guild_mf_f) : Gmf.append(b)
            # Count number of bacteria in each guild
            self._numGrf.append(len(Grf))
            self._numGrm.append(len(Grm))
            self._numGpf.append(len(Gpf))
            self._numGpm.append(len(Gpm))
            self._numGmm.append(len(Gmm))
            self._numGmf.append(len(Gmf))
            self._totalMicrobes.append(len(bacter))
            # Register the states that each guild is in.
            self._Grf_state.append(GuildState(Grf))
            self._Grm_state.append(GuildState(Grm))
            self._Gpf_state.append(GuildState(Gpf))
            self._Gpm_state.append(GuildState(Gpm))
            self._Gmm_state.append(GuildState(Gmm))
            self._Gmf_state.append(GuildState(Gmf))
            self._time.append(time)
            if nutrient_input:
                self.nutrient_inputs.append(nutrient_input)
            # if self._logCounter % (self._logPeriod * 24.0) == 0:
            #     excreted = self._vessel.bacteria_excreted_logging
            #     load = len(bacter)
            #     if load > 0:
            #         generations = float(excreted) / load
            #         print('generations = ' + str(generations))
            #     self._vessel.bacteria_excreted_logging = 0  # reset.

        # Can't do this for the very first iteration, as nothing stepped yet.
        if time >= self.bacteria_graphing_start_time and time > 0:
            # Graph_global_dynamics drawing stuff. Do this only periodically
            # ---------
            # DEPRECATED in favour of heatmap-based plots. 
            # filename = self._graphPath + '/limit-' + str(time) + '.svg'
            # self.plot_limitationCN_spatial(filename, bacter, Grf, Grm, Gpf, Gpm, Gmm, Gmf)
            # filename = self._graphPath + '/limitResource-' + str(time) + '.svg'
            # self.plot_guild_limiting_resource(filename, Grf, Grm, Gpf, Gpm, Gmm, Gmf)
            # filename = self._graphPath + '/rates-' + str(time) + '.svg'
            # self.plot_guild_rates(filename, Grf, Grm, Gpf, Gpm, Gmm, Gmf)
            # ---------
            self.heatmap_logger.log(time, self._vessel._nutrientStores, bacter)
        self._logCounter += 1

    def string_state(self):
        return ('mouse ' + str(self._vessel.get_id())  +
                '; Grf=' + str(self._numGrf[-1]) + ', Grm=' + str(self._numGrm[-1]) + ', Gpf=' + str(self._numGpf[-1]) +
                ', Gpm=' + str(self._numGpm[-1]) + ', Gmm=' + str(self._numGmm[-1]) + ', Gmf=' + str(self._numGmf[-1]) +
                '; nns=' + str(len(self._vessel._nutrientStores)))

    def plot_spatial_nutrient_distribution(self):
        """ Plots nutrient concentrations along the length of the GIT. This is only for the last time point """
        nut = self._vessel._nutrientStores
        xs = list(range(len(nut)))
        plt.clf()

        Nfs = [i.Nf for i in nut]
        Nms = [i.Nm for i in nut]
        Cws = [i.Cw for i in nut]
        Cds = [i.Cd for i in nut]
        Cms = [i.Cm for i in nut]
        Cis = [i.Ci for i in nut]
        data_df = pandas.DataFrame(  # Used to write raw data to the file system. 
            {'Location': xs,
             'Nf': Nfs,
             'Nm': Nms,
             'Cw': Cws,
             'Cd': Cds, 
             'Cm': Cms, 
             'Ci': Cis
            })
        
        p1, = plt.plot(xs, Nfs, 'r-')
        p2, = plt.plot(xs, Nms, 'g-')
        p3, = plt.plot(xs, Cws, 'm-')
        p4, = plt.plot(xs, Cds, 'c-')
        p5, = plt.plot(xs, Cms, 'b-')
        p6, = plt.plot(xs, Cis, 'y-')
        try:
            plt.yscale('log')
        except ValueError:
            print('Warning! No data to plot, likely all nutrients depleted along the GIT. Cannot create log scale.')
            plt.yscale('linear')
        plt.xlabel('Gastrointestinal tract')
        plt.ylabel('Grams')
        plt.legend(['Nf','Nm','Cw','Cd','Cm','Ci'], loc='upper center', bbox_to_anchor=(0.5, 1.1),
                   ncol=6,prop={'size':8})
        plt.rcParams.update({'font.size': 18})
        plt.savefig(self._graphPath + '/spatial_nutrients.png' ,dpi=300)
        # Write raw data to the file system
        data_df.to_csv(self._graphPath + '/spatial_nutrients.csv', index=False)

    def plot_spatial_distribution(self, filename):
        """ Draw a stacked bar chart representing guilds as a proportion of all bacteria, along length of colon """
        sbg = StackedBarGrapher()
        bacter = self._vessel._bacteria
        segments = 50  # How many bars to draw, how many segments gut length is divided into
        segWidth = int(len(bacter) / segments)  # Rounds to nearest int.
        count = segWidth  # Count down how many bacteria in the current segment have been processed.
        Grf = 0  # Count how many of each type in each segment.
        Grm = 0
        Gpf = 0
        Gpm = 0
        Gmm = 0
        Gmf = 0
        # Store the proportions of each guild in each segment. Items in list represent segments.
        propGrf = []; propGrm = []; propGpf = []; propGpm = []; propGmm = []; propGmf = []
        for b in bacter:
            # Iterate through each bacterium, log number encountered.
            if   isinstance(b, gutsim.bacteria.Guild_rf_f) and not b.is_dead(): Grf += 1
            elif isinstance(b, gutsim.bacteria.Guild_rm_m) and not b.is_dead(): Grm += 1
            elif isinstance(b, gutsim.bacteria.Guild_pf_f) and not b.is_dead(): Gpf += 1
            elif isinstance(b, gutsim.bacteria.Guild_pm_m) and not b.is_dead(): Gpm += 1
            elif isinstance(b, gutsim.bacteria.Guild_m_m) and not b.is_dead(): Gmm += 1
            elif isinstance(b, gutsim.bacteria.Guild_mf_f) and not b.is_dead(): Gmf += 1
            count -= 1
            if count == 0:  # Reached end of segment. Tally up.
                total = Grf + Grm + Gpf + Gpm + Gmm + Gmf
                if total > 0:  # In case all these bugs are dead. Hence total alive = 0.
                    propGrf.append(float(Grf) / total)
                    propGrm.append(float(Grm) / total)
                    propGpf.append(float(Gpf) / total)
                    propGpm.append(float(Gpm) / total)
                    propGmm.append(float(Gmm) / total)
                    propGmf.append(float(Gmf) / total)
                    Grf = 0; Grm = 0; Gpf = 0; Gpm = 0; Gmm = 0; Gmf = 0
                    count = segWidth

        xLabels = [''] * segments
        xLabels[0] = 'proximal'
        xLabels[-5] = 'distal'
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        # Same order for data and colours.
        d_colors = [
            guild_colours['Gmm'],
            guild_colours['Gmf'],
            guild_colours['Gpm'],
            guild_colours['Gpf'],
            guild_colours['Grm'],
            guild_colours['Grf']
        ]
        data = np.array([propGmm, propGmf, propGpm, propGpf, propGrm, propGrf])
        data = np.transpose(data)
        sbg.stackedBarPlot(ax=ax1, data=data, cols=d_colors,)
        # plt.title('Community composition along GIT')
        plt.ylabel('Proportion')
        plt.yticks([1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
        plt.xlabel('Gastrointestinal tract')
        plt.xticks(range(segments), xLabels)
        plt.legend(['Gmm', 'Gmf', 'Gpm', 'Gpf', ' Grm', ' Grf'], loc='upper center', bbox_to_anchor=(0.5, 1.1),
                    ncol=5, prop={'size' : 14})
        plt.tight_layout()
        plt.rcParams.update({'font.size': 18})
        plt.savefig(filename + ".svg", dpi=300)
        plt.close(fig)
        del fig

    def plot_spatial_distribution_bacterial_states(self, filename):
        sbg = StackedBarGrapher()
        bacter = self._vessel._bacteria
        segments = 50
        segWidth = int(len(bacter) / segments)  # Rounds to nearest int.
        count = segWidth
        exponential = 0; stressed = 0; resistant = 0; dead = 0
        # Store the proportions of each guild in each segment. Items in list represent segments.
        propExp = []; propStress = []; propRes = []; propDead = []
        for b in bacter:
            # Iterate through each bacterium, log number encountered.
            if   b.is_dead()		:   dead += 1
            elif b.is_resistant()   :   resistant += 1
            elif b.is_stressed()	:   stressed += 1
            elif b.is_exponential() :   exponential += 1

            count -= 1
            if count == 0:
                total = exponential + stressed + resistant + dead
                propExp.append(float(exponential) / total)
                propStress.append(float(stressed) / total)
                propRes.append(float(resistant) / total)
                propDead.append(float(dead) / total)
                exponential = 0; stressed = 0; resistant = 0; dead = 0
                count = segWidth

        xLabels = [''] * segments
        xLabels[0] = 'cecum'
        xLabels[-5] = 'distal colon'
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        data = np.array([propExp, propStress, propRes, propDead])
        data = np.transpose(data)
        d_colors = ['g', 'b', 'r', 'k']
        sbg.stackedBarPlot(ax=ax1, data=data, cols=d_colors,)
        #plt.title('Community composition along GIT')
        plt.ylabel('Proportion')
        plt.xlabel('Gastrointestinal tract')
        plt.xticks(range(segments), xLabels)
        plt.legend(['exp', 'stress', 'res', 'dead'],
                   loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=5, prop={'size': 14})
        plt.tight_layout()
        plt.rcParams.update({'font.size': 18})
        plt.savefig(filename + ".svg", dpi=300)
        plt.close(fig)

    def plot_microbial_loads(self, filename, title=None):
        """
        Draws a graph_global_dynamics of this mouse simulation.
        Graph shows both bacterial abundances and nutrient quantities over time.
        Graph is saved to the supplied `filename`.
        If supplied, the text in string `title` will be printed as the graph_global_dynamics title.
        """
        print(self.string_state())
        plt.clf()
        if title is not None:
            plt.title(title)
        day_time = [t / 24.0 for t in self._time ]
        
        df = pandas.DataFrame.from_dict({
            'Time (days)': day_time,
            'Grf': self._numGrf,
            'Grm': self._numGrm,
            'Gpf': self._numGpf,
            'Gpm': self._numGpm,
            'Gmf': self._numGmf,
            'Gmm': self._numGmm,
            'Total': self._totalMicrobes
        })
        # Turn into long-form dataframe for seaborn. 
        df = pandas.melt(
            df,
            id_vars=['Time (days)'],
            value_vars=['Grf', 'Grm', 'Gpf', 'Gpm', 'Gmf', 'Gmm', 'Total'],
            var_name='Guild', value_name='Microbial load')
        ax = sns.lineplot(
            x='Time (days)', y='Microbial load', hue='Guild', data=df,
            palette=guild_colours,  # palette="hls",
            legend=False)  # May wish to turn this on or off. For preparing papers it can be annoying.
        plt.savefig(filename + ".png", dpi=300)
        plt.savefig(filename + ".svg")
        
    def graph_nutrient_inputs(self, filename, title=None):
        """
        Plots the quantity of nutrient input every timestep (not rate, just quantity) as time series data.
        """
        plt.clf()
        if title is not None:
            plt.title(title)
        # Prepare data
        dayTime = [t / 24.0 for t in self._time[1:]]  # Drop first element, as time zero no nutrient have been input.
        Cw = [n.Cw * 1000 for n in self.nutrient_inputs]
        Cd = [n.Cd * 1000 for n in self.nutrient_inputs]
        Cm = [n.Cm * 1000 for n in self.nutrient_inputs]
        Ci = [n.Ci * 1000 for n in self.nutrient_inputs]
        Nf = [n.Nf * 1000 for n in self.nutrient_inputs]
        Nm = [n.Nm * 1000 for n in self.nutrient_inputs]
        # Do the plotting.
        p1, = plt.plot(dayTime, Cw, 'm-')
        p2, = plt.plot(dayTime, Cd, 'c-')
        p3, = plt.plot(dayTime, Cm, 'b-')
        p4, = plt.plot(dayTime, Nf, 'r-')
        p5, = plt.plot(dayTime, Nm, 'g-')
        p6, = plt.plot(dayTime, Ci, 'y-')
        plt.xlabel('Time (days)')
        if max(dayTime) >= 5.0:
            majorLocator   = MultipleLocator(5)  # Set major tick marks every 5 days.
            majorFormatter = FormatStrFormatter('%d')
            minorLocator   = MultipleLocator(1)  # Set minor tick marks every day
            ax = plt.gca()
            ax.xaxis.set_major_locator(majorLocator)
            ax.xaxis.set_major_formatter(majorFormatter)
            # For the minor ticks, use no labels; default NullFormatter
            ax.xaxis.set_minor_locator(minorLocator)
        plt.ylabel('Nutrient (grams)')
        plt.legend([p1, p2, p3, p4, p5, p6], ['Cw', 'Cd', 'Cm', 'Nf', 'Nm', 'Ci'], loc=2)
        plt.rcParams.update({'font.size': 18})
        plt.savefig(filename + ".svg", dpi=300)

    def graph_guild_states(self, prefix):
        """
        Plots time series stacked area graphs of the proportions of each guild that reside in each possible state.
        """
        def compile_state_data(guild):
            """ `Guild` is a list of GuildState objects, this extracts the data into a format for graphing. """
            dead = [s.stack_dead() for s in guild]
            resistant = [s.stack_resistant() for s in guild]
            stressed = [s.stack_stressed() for s in guild]
            exponential = [s.stack_exponential() for s in guild]
            return dead, resistant, stressed, exponential

        def graph_guild(filename, title, states):
            (dead, resistant, stressed, exponential) = states
            plt.clf()
            plt.title(title)
            plt.fill_between(self._time,0,dead, facecolor='black', interpolate=True)
            plt.fill_between(self._time,dead,resistant, facecolor='red', interpolate=True)
            plt.fill_between(self._time,resistant,stressed, facecolor='blue', interpolate=True)
            plt.fill_between(self._time,stressed,exponential, facecolor='green', interpolate=True)
            plt.xlabel('Time (hours)')
            plt.ylabel('Proportion of population (%)')
            plt.ylim((0,100))
            plt.savefig(filename + ".svg", dpi=300)

        def graph_guild_limit_nutrient(filename, title, guild):
            c_lim = [s.stack_c_limited() for s in guild]
            n_lim = [s.stack_n_limited() for s in guild]
            plt.clf()
            plt.title('Limiting resource')
            plt.fill_between(self._time, 0, c_lim, facecolor='red', interpolate=True)
            plt.fill_between(self._time, c_lim, n_lim, facecolor='green', interpolate=True)
            plt.xlabel('Time (hours)')
            plt.ylabel('Proportion of population (%)')
            plt.ylim((0,100))
            plt.savefig(filename + ".svg", dpi=300)

        def graph_guild_doubling_times(filename, data, guild_label, color, ylim, title=None):
            # Filter out the particular guild of interest. 
            df = data.loc[data['Guild'] == guild_label, :]
            plt.clf()  
            ax = sns.lineplot(x='Time (days)', y='Doubling time (hrs)', data=df, color=color)
            if title:
                plt.title(title)
            plt.savefig(filename + ".svg")
            plt.savefig(filename + ".png", dpi=300)


        graph_guild(prefix + '_states_Grf.eps', 'Grf guild states', compile_state_data(self._Grf_state))
        graph_guild(prefix + '_states_Grm.eps', 'Grm guild states', compile_state_data(self._Grm_state))
        graph_guild(prefix + '_states_Gpf.eps', 'Gpf guild states', compile_state_data(self._Gpf_state))
        graph_guild(prefix + '_states_Gpm.eps', 'Gpm guild states', compile_state_data(self._Gpm_state))
        graph_guild(prefix + '_states_Gmm.eps', 'Gmm guild states', compile_state_data(self._Gmm_state))
        graph_guild(prefix + '_states_Gmf.eps', 'Gmf guild states', compile_state_data(self._Gmf_state))
        # Plot area graph of which nutrient is limiting across each guild over time.
        sns.set(style='white', context='talk')
        graph_guild_limit_nutrient(prefix + '_limNut_Grf', 'Grf limiting resources', self._Grf_state)
        graph_guild_limit_nutrient(prefix + '_limNut_Grm', 'Grm limiting resources', self._Grm_state)
        graph_guild_limit_nutrient(prefix + '_limNut_Gpf', 'Gpf limiting resources', self._Gpf_state)
        graph_guild_limit_nutrient(prefix + '_limNut_Gpm', 'Gpm limiting resources', self._Gpm_state)
        graph_guild_limit_nutrient(prefix + '_limNut_Gmm', 'Gmm limiting resources', self._Gmm_state)
        graph_guild_limit_nutrient(prefix + '_limNut_Gmf', 'Gmf limiting resources', self._Gmf_state)
        sns.set(style='darkgrid', context='talk')
        # Plot the doubling time for the guilds.
        plt.clf()
        time_col = []
        doubling_time_col = []
        guild_col = []
        guild_state_collection = [self._Grf_state, self._Grm_state, 
                                 self._Gpf_state, self._Gpm_state, 
                                 self._Gmm_state, self._Gmf_state]
        guild_labels = ['Grf', 'Grm', 'Gpf', 'Gpm', 'Gmm', 'Gmf']
        guild_colours = sns.color_palette("hls", n_colors=len(guild_labels))
        # Compile each column of the long-form DataFrame. 
        for guild_states, guild_label in zip(guild_state_collection, guild_labels):
            for t, state in zip(self._time, guild_states):
                n_observations = len(state.doubling_times) 
                t_days = t / 24.                   
                time_col.extend([t_days] * n_observations)
                doubling_time_col.extend(state.doubling_times)
                guild_col.extend([guild_label] * n_observations)

        df = pandas.DataFrame.from_dict({'Time (days)': time_col, 'Guild': guild_col, 
                                         'Doubling time (hrs)': doubling_time_col})
        ax = sns.lineplot(x='Time (days)', y='Doubling time (hrs)', hue='Guild', data=df)        
        plt.savefig(prefix + '_doubling_times.svg')
        plt.savefig(prefix + '_doubling_times.png', dpi=300)

        # Set common y axis labels for all guilds. 
        # all_states = self._Grf_state[1:] + self._Grm_state[1:] + self._Gpf_state[1:] + self._Gpm_state[1:] + self._Gmm_state[1:] + self._Gmf_state[1:]       
        # max_doubling_time = np.max([np.max(gs.doubling_times) for gs in all_states if len(gs.doubling_times) > 0])
        y_axis_lims = [0, df.loc[:, 'Doubling time (hrs)'].max()]
        graph_guild_doubling_times(filename=prefix + '_doubling_times_Grf', data=df, guild_label='Grf',
                                   color=guild_colours[0], ylim=y_axis_lims)
        graph_guild_doubling_times(filename=prefix + '_doubling_times_Grm', data=df, guild_label='Grm',
                                   color=guild_colours[1], ylim=y_axis_lims)
        graph_guild_doubling_times(filename=prefix + '_doubling_times_Gpf', data=df, guild_label='Gpf',
                                   color=guild_colours[2], ylim=y_axis_lims)
        graph_guild_doubling_times(filename=prefix + '_doubling_times_Gpm', data=df, guild_label='Gpm',
                                   color=guild_colours[3], ylim=y_axis_lims)
        graph_guild_doubling_times(filename=prefix + '_doubling_times_Gmm', data=df, guild_label='Gmm',
                                   color=guild_colours[4], ylim=y_axis_lims)
        graph_guild_doubling_times(filename=prefix + '_doubling_times_Gmf', data=df, guild_label='Gmf',
                                   color=guild_colours[5], ylim=y_axis_lims)            

    def relativeAbundances(self, index=-1):
        """
        Calculates the relative abundances of bacteria logged, based on the last data submitted to the logger, or a
        specified index (time point) if specified.
        Relative abundances are returned as a dictionary, indexed by the bacteria code.
        """
        # Last time point unless other index is specified.
        total = float(self._numGrf[index] + self._numGrm[index] + self._numGpf[index] + self._numGpm[index] +
                      self._numGmm[index] + self._numGmf[index])
        abundances = {}
        if total != 0:
            abundances['Grf'] = (100.0 * float(self._numGrf[index])) / total
            abundances['Grm'] = (100.0 * float(self._numGrm[index])) / total
            abundances['Gpf'] = (100.0 * float(self._numGpf[index])) / total
            abundances['Gpm'] = (100.0 * float(self._numGpm[index])) / total
            abundances['Gmm'] = (100.0 * float(self._numGmm[index])) / total
            abundances['Gmf'] = (100.0 * float(self._numGmf[index])) / total
        else:  # Safety, avoid div by zero.
            abundances['Grf'] = 0.0
            abundances['Grm'] = 0.0
            abundances['Gpf'] = 0.0
            abundances['Gpm'] = 0.0
            abundances['Gmm'] = 0.0
            abundances['Gmf'] = 0.0
        return abundances

    def writeMouseTimeseriesToFile(self, filename):
        """ Will write Mouse timeseries data to the specified file. """
        f = open(filename,'w')

        f.write('#{:8s},{:16s},{:16s},'
                '{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},'
                '{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},'
                '{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},'
                '{:9s},'
                '{:7s},{:7s},{:7s},{:7s}\n'
            .format('time', 'carb-intake-KJ/d', 'prot-intake-KJ/d',
                'relGrf', 'relGrm', 'relGpf', 'relGpm', 'relGmm','relGmf',
                'absGrf', 'absGrm', 'absGpf', 'absGpm', 'absGmm', 'absGmf', 'total',
                'nLimGrf', 'nLimGrm', 'nLimGpf', 'nLimGpm', 'nLimGmm', 'nLimGmf', 'extime',
                'dead', 'resist', 'stress', 'expon'))
        ci = self._vessel.daily_wheatstarch_kj + self._vessel.daily_dextrinised_kj + self._vessel.daily_sucrose_kj
        pi = self._vessel.daily_nitrogen_kj
        for i, t in enumerate(self._time):
            relab = self.relativeAbundances(index=i)
            relGrf = relab['Grf']
            relGrm = relab['Grm']
            relGpf = relab['Gpf']
            relGpm = relab['Gpm']
            relGmm = relab['Gmm']
            relGmf = relab['Gmf']
            Grf = self._numGrf[i]
            Grm = self._numGrm[i]
            Gpf = self._numGpf[i]
            Gpm = self._numGpm[i]
            Gmm = self._numGmm[i]
            Gmf = self._numGmf[i]
            GrfNlim = self._Grf_state[i].n_limited
            GrmNlim = self._Grm_state[i].n_limited
            GpfNlim = self._Gpf_state[i].n_limited
            GpmNlim = self._Gpm_state[i].n_limited
            GmmNlim = self._Gmm_state[i].n_limited
            GmfNlim = self._Gmf_state[i].n_limited
            total = Grf + Grm + Gpf + Gpm + Gmm + Gmf
            execution_time = (time.time() - self._start_time) / 60.0  # Convert seconds to min.
            all_dead = [b for b in self._vessel._bacteria if b.is_dead()]
            all_resistant = [b for b in self._vessel._bacteria if b.is_resistant()]
            all_stressed = [b for b in self._vessel._bacteria if b.is_stressed()]
            all_exponential = [b for b in self._vessel._bacteria if b.is_exponential()]
            # Relative abundance of bacteria in each state amongst all bacteria
            all_dead_prop = 100.0 * float(len(all_dead)) / float(total)
            all_resistant_prop = 100.0 * float(len(all_resistant)) / float(total)
            all_stressed_prop = 100.0 * float(len(all_stressed)) / float(total)
            all_exponential_prop = 100.0 * float(len(all_exponential)) / float(total)
            f.write('{:9f},{:16.4f},{:16.4f},'   # vessel ID, ci and pi
                    '{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},'   # relative abundances
                    '{:7d},{:7d},{:7d},{:7d},{:7d},{:7d},{:7d},'   # abundances, total bacteria load
                    '{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},'   # proportion of each guild N limited
                    '{:9.0f},'   # execution time
                    '{:7.2f},{:7.2f},{:7.2f},{:7.2f}\n'   # bacteria states
                .format(t, ci, pi,
                            relGrf, relGrm, relGpf, relGpm, relGmm, relGmf,
                            Grf, Grm, Gpf, Gpm, Gmm, Gmf, total,
                            GrfNlim, GrmNlim, GpfNlim, GpmNlim, GmmNlim, GmfNlim,
                            execution_time,
                            all_dead_prop, all_resistant_prop, all_stressed_prop, all_exponential_prop))
        f.close()

    def writeMouseToFile(self, filename, xml_filename):
        """
        Will write Mouse data to the specified file.
        """
        def log_output(parent, name, data):
            child = SubElement(parent, name)
            child.text = str(data)

        def prettify(elem):
            """
            Return a pretty-printed XML string for the Element.

            Taken from: http://pymotw.com/2/xml/etree/ElementTree/create.html on 19/08/2014
            """
            rough_string = ElementTree.tostring(elem, 'utf-8')
            reparsed = minidom.parseString(rough_string)
            return reparsed.toprettyxml(indent="  ")

        f = open(filename,'w')

        f.write('#{:8s},{:16s},{:16s},'
                '{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},'
                '{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},'
                '{:7s},{:7s},{:7s},{:7s},{:7s},{:7s},'
                '{:9s},'
                '{:7s},{:7s},{:7s},{:7s}\n'
            .format('mouseID', 'carb-intake-KJ/d', 'prot-intake-KJ/d',
                'relGrf', 'relGrm', 'relGpf', 'relGpm', 'relGmm','relGmf',
                'absGrf', 'absGrm', 'absGpf', 'absGpm', 'absGmm', 'absGmf', 'total',
                'nLimGrf', 'nLimGrm', 'nLimGpf', 'nLimGpm', 'nLimGmm', 'nLimGmf', 'extime',
                'dead', 'resist', 'stress', 'expon'))
        ci = self._vessel.daily_wheatstarch_kj + self._vessel.daily_dextrinised_kj + self._vessel.daily_sucrose_kj
        pi = self._vessel.daily_nitrogen_kj
        relab = self.relativeAbundances()
        relGrf = relab['Grf']
        relGrm = relab['Grm']
        relGpf = relab['Gpf']
        relGpm = relab['Gpm']
        relGmm = relab['Gmm']
        relGmf = relab['Gmf']
        Grf = self._numGrf[-1]
        Grm = self._numGrm[-1]
        Gpf = self._numGpf[-1]
        Gpm = self._numGpm[-1]
        Gmm = self._numGmm[-1]
        Gmf = self._numGmf[-1]
        GrfNlim = self._Grf_state[-1].n_limited
        GrmNlim = self._Grm_state[-1].n_limited
        GpfNlim = self._Gpf_state[-1].n_limited
        GpmNlim = self._Gpm_state[-1].n_limited
        GmmNlim = self._Gmm_state[-1].n_limited
        GmfNlim = self._Gmf_state[-1].n_limited
        total = Grf + Grm + Gpf + Gpm + Gmm + Gmf
        execution_time = (time.time() - self._start_time) / 60.0   # convert seconds to min.
        all_dead = [b for b in self._vessel._bacteria if b.is_dead()]
        all_resistant = [b for b in self._vessel._bacteria if b.is_resistant()]
        all_stressed = [b for b in self._vessel._bacteria if b.is_stressed()]
        all_exponential = [b for b in self._vessel._bacteria if b.is_exponential()]
        # relative abundance of bacteria in each state amongst all bacteria
        all_dead_prop = 100.0 * float(len(all_dead)) / float(total)
        all_resistant_prop = 100.0 * float(len(all_resistant)) / float(total)
        all_stressed_prop = 100.0 * float(len(all_stressed)) / float(total)
        all_exponential_prop = 100.0 * float(len(all_exponential)) / float(total)
        f.write('{:9d},{:16.4f},{:16.4f},'   # vessel ID, ci and pi
                '{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},'   # relative abundances
                '{:7d},{:7d},{:7d},{:7d},{:7d},{:7d},{:7d},'   # abundances, total bacteria load
                '{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},{:7.2f},'   # proportion of each guild N limited
                '{:9.0f},'   # execution time
                '{:7.2f},{:7.2f},{:7.2f},{:7.2f}\n'   # bacteria states
                .format(self._vessel._ID, ci, pi,
                        relGrf, relGrm, relGpf, relGpm, relGmm, relGmf,
                        Grf, Grm, Gpf, Gpm, Gmm, Gmf, total,
                        GrfNlim, GrmNlim, GpfNlim, GpmNlim, GmmNlim, GmfNlim,
                        execution_time,
                        all_dead_prop, all_resistant_prop, all_stressed_prop, all_exponential_prop))
        f.close()
        if False:
            root = Element('output')
            log_output(root, 'mouseID', self._vessel._ID)
            log_output(root, 'carb-intake-KJ_d', ci)
            log_output(root, 'prot-intake-KJ_d', pi)
            log_output(root, 'relGrf', relGrf)
            log_output(root, 'relGrm', relGrm)
            log_output(root, 'relGpf', relGpf)
            log_output(root, 'relGpm', relGpm)
            log_output(root, 'relGmm', relGmm)
            log_output(root, 'absGrf', Grf)
            log_output(root, 'absGrm', Grm)
            log_output(root, 'absGpf', Gpf)
            log_output(root, 'absGpm', Gpm)
            log_output(root, 'absGmm', Gmm)
            log_output(root, 'total', total)
            log_output(root, 'nLimGrf', GrfNlim)
            log_output(root, 'nLimGrm', GrmNlim)
            log_output(root, 'nLimGpf', GpfNlim)
            log_output(root, 'nLimGpm', GpmNlim)
            log_output(root, 'nLimGmm', GmmNlim)
            with open(xml_filename, 'w') as f:
                f.write(prettify(root))

    def plot_spatiotemporal_heatmaps(self):
        self.heatmap_logger.write(self._graphPath)

    def plotBacterialNutrients(self, filename, bacter=[], Grf=[], Grm=[], Gpf=[], Gpm=[], Gmm=[], Gmf=[]):
        GrfC = []; GrmC = []; GpfC = []; GpmC = []; GmmC = []; GmfC = []
        GrfN = []; GrmN = []; GpfN = []; GpmN = []; GmmN = []; GmfN = []
        if Grf:
            for b in Grf:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                GrfC.append(carbon)
                GrfN.append(nitrogen)
            for b in Grm:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                GrmC.append(carbon)
                GrmN.append(nitrogen)
            for b in Gpf:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                GpfC.append(carbon)
                GpfN.append(nitrogen)
            for b in Gpm:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                GpmC.append(carbon)
                GpmN.append(nitrogen)
            for b in Gmm:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                GmmC.append(carbon)
                GmmN.append(nitrogen)
            for b in Gmf:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                GmfC.append(carbon)
                GmfN.append(nitrogen)
        else:       # in case bacteria have not already been separated.
            for b in bacter:
                (carbon, nitrogen) = b.get_carbon_nitrogen()
                if isinstance(b, gutsim.bacteria.Guild_rf_f):      GrfC.append(carbon), GrfN.append(nitrogen)
                elif isinstance(b, gutsim.bacteria.Guild_rm_m):    GrmC.append(carbon), GrmN.append(nitrogen)
                elif isinstance(b, gutsim.bacteria.Guild_pf_f):    GpfC.append(carbon), GpfN.append(nitrogen)
                elif isinstance(b, gutsim.bacteria.Guild_pm_m):    GpmC.append(carbon), GpmN.append(nitrogen)
                elif isinstance(b, gutsim.bacteria.Guild_m_m):    GmmC.append(carbon), GmmN.append(nitrogen)
                elif isinstance(b, gutsim.bacteria.Guild_mf_f):    GmfC.append(carbon), GmfN.append(nitrogen)
        plt.clf()
        plt.hold(True)

        # plot the data in pairs, and change the colors of boxes as you go.
        bp = plt.boxplot([GrfC, GrfN], positions=[1,2], widths=0.5)
        setBoxColors(bp)
        bp = plt.boxplot([GrmC, GrmN], positions=[3,4], widths=0.5)
        setBoxColors(bp)
        bp = plt.boxplot([GpfC, GpfN], positions=[5,6], widths=0.5)
        setBoxColors(bp)
        bp = plt.boxplot([GpmC, GpmN], positions=[7,8], widths=0.5)
        setBoxColors(bp)
        bp = plt.boxplot([GmmC, GmmN], positions=[9,10], widths=0.5)
        setBoxColors(bp)
        bp = plt.boxplot([GmfC, GmfN], positions=[11,12], widths=0.5)
        setBoxColors(bp)

        #plt.yscale('log')
        plt.xlim([0,13])
        #plt.ylim([0.1,100])
        plt.xticks(range(1,10+1), ['GrfC', 'GrfN', 'GrmC', 'GrmN', 'GpfC', 'GpfN', 'GpmC', 'GpmN', 'GmmC', 'GmmN', 'GmfC', 'GmfN'])
        plt.xlabel('Guild & nutrient')
        plt.ylabel('Quantity')
        plt.savefig(filename)


    def plot_limitationCN_spatial(self, filename, bacter=[], Grf=[], Grm=[], Gpf=[], Gpm=[], Gmm=[], Gmf=[]):
        """
        Creates a stacked bar chart, representing the proportion of bacteria in each guild that are carbon limited,
        or nitrogen limited.
        """
        # will store how many bacteria there are of each guild which are either limited by carbon or nitrogen in these.
        GrfC = 0.0; GrmC = 0.0; GpfC = 0.0; GpmC = 0.0; GmmC = 0.0; GmfC = 0.0
        GrfN = 0.0; GrmN = 0.0; GpfN = 0.0; GpmN = 0.0; GmmN = 0.0; GmfN = 0.0
        if Grf:
            for b in Grf:
                if b.carbon_limiting():		 GrfC += 1.0
                elif b.nitrogen_limiting():	 GrfN += 1.0
            for b in Grm:
                if b.carbon_limiting():		 GrmC += 1.0
                elif b.nitrogen_limiting():	 GrmN += 1.0
            for b in Gpf:
                if b.carbon_limiting():		 GpfC += 1.0
                elif b.nitrogen_limiting():	 GpfN += 1.0
            for b in Gpm:
                if b.carbon_limiting():		 GpmC += 1.0
                elif b.nitrogen_limiting():	 GpmN += 1.0
            for b in Gmm:
                if b.carbon_limiting():		 GmmC += 1.0
                elif b.nitrogen_limiting():	 GmmN += 1.0
            for b in Gmf:
                if b.carbon_limiting():		 GmfC += 1.0
                elif b.nitrogen_limiting():	 GmfN += 1.0
        else:
            for b in bacter:
                if isinstance(b, gutsim.bacteria.Guild_rf_f):
                    if b.carbon_limiting():		 GrfC += 1.0
                    elif b.nitrogen_limiting():	 GrfN += 1.0
                elif isinstance(b, gutsim.bacteria.Guild_rm_m):
                    if b.carbon_limiting():		 GrmC += 1.0
                    elif b.nitrogen_limiting():	 GrmN += 1.0
                elif isinstance(b, gutsim.bacteria.Guild_pf_f):
                    if b.carbon_limiting():		 GpfC += 1.0
                    elif b.nitrogen_limiting():	 GpfN += 1.0
                elif isinstance(b, gutsim.bacteria.Guild_pm_m):
                    if b.carbon_limiting():		 GpmC += 1.0
                    elif b.nitrogen_limiting():	 GpmN += 1.0
                elif isinstance(b, gutsim.bacteria.Guild_m_m):
                    if b.carbon_limiting():		 GmmC += 1.0
                    elif b.nitrogen_limiting():	 GmmN += 1.0
                elif isinstance(b, gutsim.bacteria.Guild_mf_f):
                    if b.carbon_limiting():		 GmfC += 1.0
                    elif b.nitrogen_limiting():	 GmfN += 1.0
        # Convert the raw data into percentages. Some type safety here in case all of a bacteria species goes extinct.
        total = GrfC + GrfN
        if total != 0.0:
            GrfC = (GrfC / total) * 100.0
            GrfN = (GrfN / total) * 100.0
        else: GrfC = GrfN = 0.0

        total = GrmC + GrmN
        if total != 0.0:
            GrmC = (GrmC / total) * 100.0
            GrmN = (GrmN / total) * 100.0
        else: GrmC = GrmN = 0.0

        total = GpfC + GpfN
        if total != 0.0:
            GpfC = (GpfC / total) * 100.0
            GpfN = (GpfN / total) * 100.0
        else: GpfC = GpfN = 0.0

        total = GpmC + GpmN
        if total != 0.0:
            GpmC = (GpmC / total) * 100.0
            GpmN = (GpmN / total) * 100.0
        else: GpmC = GpmN = 0.0

        total = GmmC + GmmN
        if total != 0.0:
            GmmC = (GmmC / total) * 100.0
            GmmN = (GmmN / total) * 100.0
        else: GmmC = GmmN = 0.0

        total = GmfC + GmfN
        if total != 0.0:
            GmfC = (GmfC / total) * 100.0
            GmfN = (GmfN / total) * 100.0
        else: GmfC = GmfN = 0.0

        ind = np.arange(6)
        width = 0.7	   # the width of the bars
        C = [GrfC, GrmC, GpfC, GpmC, GmmC, GmfC]
        N = [GrfN, GrmN, GpfN, GpmN, GmmN, GmfN]
        plt.clf()
        p1 = plt.bar(ind, C, width, color='r')
        p2 = plt.bar(ind, N, width, color='g', bottom=C)
        plt.ylim([0,100])
        plt.ylabel('Percentage')
        plt.xlabel('Guild')
        plt.xticks(ind + width / 2.0, ('Grf', 'Grm', 'Gpf', 'Gpm', 'Gmm', 'Gmf'))
        plt.legend((p1[0], p2[0]), ('Carbon', 'Nitrogen'))
        plt.savefig(filename)

    def plot_guild_limiting_resource(self, filename, Grf=[], Grm=[], Gpf=[], Gpm=[], Gmm=[], Gmf=[]):
        """
        Plots the distribution of internal resources, expressed in terms of the limiting resource for each bacterial
        guild
        """
        GrfL = [b._calculate_limiting_resource() for b in Grf]
        GrmL = [b._calculate_limiting_resource() for b in Grm]
        GpfL = [b._calculate_limiting_resource() for b in Gpf]
        GpmL = [b._calculate_limiting_resource() for b in Gpm]
        GmmL = [b._calculate_limiting_resource() for b in Gmm]
        GmfL = [b._calculate_limiting_resource() for b in Gmf]
        plt.clf()
        plt.boxplot((GrfL, GrmL, GpfL, GpmL, GmmL, GmfL))
        plt.ylabel('Limiting resource')
        plt.xlabel('Guild')
        plt.xticks(np.arange(6) + 1, ('Grf', 'Grm', 'Gpf', 'Gpm', 'Gmm', 'Gmf'))
        plt.ylim((0.0, 200.0))
        #plt.yscale('log')
        plt.savefig(filename)

        # Plot using seaborn

    def plot_guild_rates(self, filename, Grf = [], Grm = [], Gpf = [], Gpm = [], Gmm = [], Gmf=[]):
        """ Plot the distribution of death and growth rates in each guild at a particular point in time. """
        # Death rates.
        GrfD = [b.death_rate for b in Grf]
        GrmD = [b.death_rate for b in Grm]
        GpfD = [b.death_rate for b in Gpf]
        GpmD = [b.death_rate for b in Gpm]
        GmmD = [b.death_rate for b in Gmm]
        GmfD = [b.death_rate for b in Gmf]
        # Division (growth) rates
        GrfG = [b.divide_rate for b in Grf]
        GrmG = [b.divide_rate for b in Grm]
        GpfG = [b.divide_rate for b in Gpf]
        GpmG = [b.divide_rate for b in Gpm]
        GmmG = [b.divide_rate for b in Gmm]
        GmfG = [b.divide_rate for b in Gmf]
        plt.clf()
        # Plot the data in pairs, and change the colors of boxes as you go.
        bp = plt.boxplot([GrfG, GrfD], positions=[1, 2], widths=0.6)
        setBoxColors(bp)
        bp = plt.boxplot([GrmG, GrmD], positions=[3, 4], widths=0.6)
        setBoxColors(bp)
        bp = plt.boxplot([GpfG, GpfD], positions=[5, 6], widths=0.6)
        setBoxColors(bp)
        bp = plt.boxplot([GpmG, GpmD], positions=[7, 8], widths=0.6)
        setBoxColors(bp)
        bp = plt.boxplot([GmmG, GmmD], positions=[9, 10], widths=0.6)
        setBoxColors(bp)
        bp = plt.boxplot([GmfG, GmfD], positions=[11, 12], widths=0.6)
        setBoxColors(bp)
        # annotate the graphs
        plt.yscale('log')
        plt.ylabel('Probability')
        plt.xlabel('Guild')
        plt.ylim((0.001, 1.0))
        plt.xlim((0, 13))
        plt.xticks(range(1, 12 + 1), ['GrfG', 'GrfD',
                                      'GrmG', 'GrmD',
                                      'GpfG', 'GpfD',
                                      'GpmG', 'GpmG',
                                      'GmmG', 'GmmD',
                                      'GmfG', 'GmfD'],
                   rotation=45)
        plt.savefig(filename)



