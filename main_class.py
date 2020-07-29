######################################################################################################
######################################################################################################

# Main class for tracking and compiling all relevant data from the WRF data sets.
# In terms of tracking, this handles time, locations, and reflectivity for each cell

from netCDF4 import Dataset
from datetime import timedelta, datetime
import os
from wrf import dbz, udhel, destagger, cape_2d, srhel
from wrf import destagger, interpz3d, interplevel
import numpy as np
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature
import cartopy.crs as crs
import gc
import PyartGridder_new as PyartGridder
import pprint
from tint import Cell_tracks
from tools import hour_check, predict_location
import xarray
import pandas as pd


def get_closest_lat_lon(lat, lon, chosen_lat, chosen_lon):

        a = abs(lat-chosen_lat) + abs(lon-chosen_lon)
        i,j = np.unravel_index(a.argmin(), a.shape)

        return i,j

# Gets nearest max with in 5 grid points. Used for finding UH max near dBz max.
def get_nearby_max(uphel, i, j):
        max_uh, max_i, max_j = 0., i, j
        for x in range(-4,5):
                for y in range(-4,5):
                    try:
                            uh = uphel[i+x, j+y]

                    except IndexError:
                            continue

                    if uh > max_uh: 
                            max_uh, max_i, max_j = uh, i+x, j+y
		   
        return max_uh, max_i, max_j

# Grabs reflecitivty file. 
def get_radar_file(date, model):
        
        file_name = 'wrf2d_d01_%s_REFL_10CM_' % model
        year, month, day, hour = date.split('-')

        if int(month) == 1 or int(month) == 2 or int(month) == 3: month1, month2 = '01', '03'
        if int(month) == 4 or int(month) == 5 or int(month) == 6: month1, month2 = '04', '06'
        if int(month) == 7 or int(month) == 8 or int(month) == 9: month1, month2 = '07', '09'
        if int(month) == 10 or int(month) == 11 or int(month) == 12: month1, month2 = '10', '12'

        date_name = str(year) + str(month1) + '-' + str(year) + str(month2) + '.nc'
        return file_name + date_name
        

# Extract Class.
# This class is designed to be used to extract a single time field for whatever
# parameter you want for the cell in question.

class Extract:
    def __init__(self, filename, variable, datetime_, run_type):
        self.filename = filename
        self.variable = variable
        self.time = []
        self.time_index = ''
        self.datetime = datetime_
        self.time_i = ''
        self.run_type = run_type

        if self.run_type in self.filename:

            date_split = self.datetime.split('-')
            
            new_date = datetime(year=int(date_split[0]),month=int(date_split[1]),day=int(date_split[2]),hour=int(date_split[3]),minute=0,second=0)
            
            self.netcdf_time = np.array(xarray.open_dataset(self.filename)['Times'][...])

            self.time_i = np.where(np.array([datetime(year=int(ds.decode('utf-8').split('-')[0]),month=int(ds.decode('utf-8').split('-')[1]),day=int(ds.decode('utf-8').split('-')[2][:2]),hour=int(ds.decode('utf-8').split('-')[2][3:5]),minute=0,second=0) for ds in self.netcdf_time]) == new_date)[0][0]
            
            if 'U' in self.variable:
               
               self.cos_, self.sin_ =  xarray.open_dataset(self.filename)['COSALPHA'][...], xarray.open_dataset(self.filename)['SINALPHA'][...]
               self.lat = xarray.open_dataset(self.filename)['XLAT'][...]
               self.lon = xarray.open_dataset(self.filename)['XLONG'][...]
            
            self.param = xarray.open_dataset(self.filename)[self.variable][self.time_i,...]

        if 'constants' in self.filename:
            self.dataset = Dataset(self.filename)
            self.netcdf_time = self.dataset.variables['Time']

            self.convert_time()
            self.get_timeindex()

            self.param =  xarray.open_dataset(self.filename)[self.variable]

        self.dataset = None

    def convert_time(self):
        for hour in self.netcdf_time:
            start = datetime(year=1901,month=1,day=1,hour=0)
            delta = timedelta(hours=1) * int(hour)
            offset = start + delta
            single_time = offset.strftime('%Y-%m-%d-%H')
            parse_single_time = single_time.split('-')
            single_time = str(int(parse_single_time[0])) + '-' + str(int(parse_single_time[1])) + '-' + str(int(parse_single_time[2])) + '-' +str(int(parse_single_time[3]))

            self.time.append(single_time)

    def get_timeindex(self):
       
        for z, t in enumerate(self.time):

            if t == self.datetime: self.time_i = z

        return self.time_i
        
        
       
# This is the class that contains all of the functions for calculations involving variables from the Extract class.
# The date parser function is also here.

class Instance_Functions:
    def __init__(self, instance_list, spec_date=None):
        self.instance_list = instance_list

        # Auto run the radar and updraft helcity
        date = spec_date.split('-')
        hour = int(date[-1])

        if hour % 3 == 0:
                        

            self.converted_height = self.instance_list['Z'].param[:,...] - self.instance_list['HGT'].param[0,...]
            self.zz = xarray.DataArray(np.array(destagger(self.instance_list['Z'].param[:,...],0)) - np.array(self.instance_list['HGT'].param[0,...]))

            self.spec_date = spec_date
            self.magnitude = 0.
            self.lat = self.instance_list['U'].lat[...]
            self.lon = self.instance_list['U'].lon[...]

            self.radar = ''
            self.up_hel = ''
            self.get_date()
            self.updraft_hel_calc()
            self.bunkers_calc()


    def updraft_hel_calc(self):
        u_destag = np.array(destagger(self.instance_list['U'].param[:,:,:],2) ) 
        v_destag = np.array(destagger(self.instance_list['V'].param[:,:,:],1) )

        u_adj = u_destag * np.array(self.instance_list['U'].cos_) - v_destag * np.array(self.instance_list['U'].sin_) 
        v_adj = v_destag * np.array(self.instance_list['U'].cos_) + u_destag * np.array(self.instance_list['U'].sin_) 

        u_destag, v_destag = None, None
        self.instance_list['V'], self.instance_list['U'] = None, None

        self.u_wind = xarray.DataArray(u_adj[:,...])
        self.v_wind = xarray.DataArray(v_adj[:,...])   

        u_adj, v_adj = None, None       
        
        uphel = udhel(self.converted_height[:,:,:], 
                            self.instance_list['MAPFAC_M'].param[0,...], self.u_wind, self.v_wind, 
                            self.instance_list['W'].param[:,...], 
                            4000., 4000., bottom=2000., top=5000.)
        
        self.up_hel = np.array(uphel)
        uphel, self.instance_list['W'], self.converted_height = None, None, None



    def bunkers_calc(self):
        d = 7.5 # SharpPy default, Bunkers 2000 paper shows with smallest error. 
        
        u_wind_6km = interplevel(self.u_wind, self.zz, 6000.)
        v_wind_6km = interplevel(self.v_wind, self.zz, 6000.)

        u_wind_0km = self.u_wind[0,...]
        v_wind_0km = self.v_wind[0,...]


        mean_u_wind = mean_wind(self.zz, self.u_wind, (0,6000.))
        mean_v_wind = mean_wind(self.zz, self.v_wind, (0,6000.))

        # Bunkers calc from SharpPy.
        u_shear_6km = u_wind_6km - u_wind_0km
        v_shear_6km = v_wind_6km - v_wind_0km 

        tmp = d / (((u_shear_6km**2.) + (v_shear_6km**2.))**0.5)
        rstu = mean_u_wind + (tmp * v_shear_6km)
        rstv = mean_v_wind - (tmp * u_shear_6km)

        self.u_wind, self.v_wind = None, None
        self.zz = None

        self.mean_wind = (mean_u_wind, mean_v_wind)

        self.bunkers_uv = (rstu, rstv)
        self.instance_list['U'], self.instance_list['V'] = None, None

    def get_date(self):
        if self.spec_date:
            for key in self.instance_list:
                for i, time_item in enumerate(self.instance_list[key].time):
                    if str(self.spec_date) in str(time_item): 
                        self.instance_list[key].time_index = i
                        
        elif self.spec_date is None:
             for key in self.instance_list:
                for i, time_item in enumerate(self.instance_list[key].time):
                    self.instance_list[key].time_index = ':'

# Calculate mean wind for a layer. 
# Layer here is a tuple
def mean_wind(heights, array3d, layer):

        bottom = array3d[0,...]
        mid = interplevel(array3d, heights,(layer[1]-layer[0]) / 2.)
        top = interplevel(array3d, heights, layer[1])

        array3d_mean = (top + mid + bottom) / 3.

        return array3d_mean
            
                    
# This is the  class that puts everything together. Includes dealing with the
# awkward directory system that the NCAR WRF dataset uses.
class Instance_Dict_Builder:
    # Date needs to be input as YEAR-MM-DD-HH
    def __init__(self, target_date, run_type):
        self.target_date1 = target_date
        self.inst_list = {}
        self.run_type = run_type

        self.year = int(self.target_date1.split('-')[0])
        self.month = int(self.target_date1.split('-')[1])
        self.day = int(self.target_date1.split('-')[2])
        self.hour = int(self.target_date1.split('-')[3])
   
        self.target_date = str(self.year) + '-' + str(self.month) + '-' + str(self.day) + '-' + str(self.hour) 
        if self.month < 10:
                new_month = '0'+str(self.month)
        else:
                new_month = self.month

        if self.day < 10:
                new_day = '0'+str(self.day)
        else:
                new_day = self.day
                

        f = "%Y/%m/%d"
        d_str = str(self.year) + '/' + str(new_month) + '/' + str(new_day)
        self.compact_date = datetime.strptime(d_str, f)

        
        self.file_path = '/glade/collections/rda/data/ds612.0/%s3D/%s/' % (self.run_type, self.year)
        self.three_month_path = '/glade/collections/rda/data/ds612.0/%s/%s/' % (self.run_type, self.year)
        self.host_path = '/glade/collections/rda/data/ds612.0/INVARIANT/'
        self.radar_path = '/glade/collections/rda/data/ds612.0/%sradrefl/REFL/' % self.run_type


        self.dict_builder()

    # This is the crucial function. Calls Extract to build the variables.
    # The var_list is only needed since we caluculate updraft helicity, bunkers, mean wind
    def dict_builder(self):
        var_list = ['U','V','Z', 'W', 'P', 'Z', 'TK']
        ref_list = ['REFL_10CM']
        meta_list = ['MAPFAC_M', 'HGT']

        if float(self.target_date[-2:]) % 3 == 0:
                for thing in os.listdir(self.file_path):
                    for item in var_list:
                        if str('_' + item + '_') in str(thing):
                            file_date = thing.split('_')[4][:-3]

                            f = "%Y-%m-%d"
                            d_str = str(file_date[0:4]) + '-' + str(file_date[4:6]) + '-' +\
                                str(file_date[6:8])
                            self.file_compact_date = datetime.strptime(d_str, f) 
                                
                            if self.file_compact_date == self.compact_date:    
                                
                                data = Extract(str(self.file_path + thing), item,\
                                    self.target_date, self.run_type)          
                                self.inst_list[item] = data

        for thing in os.listdir(self.radar_path):
            for var in ref_list:
                file_date = thing.split('_')

                if '-' in file_date[4]:        

                    d = file_date[4][:-3].split('-')[0]
                    f = "%Y-%m"
                    d_str = str(d[0:4]) + '-' + str(d[4:6])
                    self.file_compact_date1 = datetime.strptime(d_str, f)

                    d = file_date[4][:-3].split('-')[1]
                    f = "%Y-%m"
                    if int(d[4:6]) == 12:
                        d_str = str(int(d[0:4])+1) + '-' + str(1)
                    else:
                        d_str = str(d[0:4]) + '-' + str(int(d[4:6])+1)
                    self.file_compact_date2 = datetime.strptime(d_str, f) 

                    
                    
                    if str('_' + var + '_') in str(thing) and \
                            self.compact_date >= self.file_compact_date1 and self.compact_date < self.file_compact_date2:  

                        data = Extract(str(self.radar_path + thing), var, self.target_date, self.run_type)          
                        self.inst_list[var] = data

      
        for thing in os.listdir(self.three_month_path):
            for var in var3_list:
                file_date = thing.split('_')

                
                if '-' in file_date[4]:        

                    d = file_date[4][:-3].split('-')[0]
                    f = "%Y-%m"
                    d_str = str(d[0:4]) + '-' + str(d[4:6])
                    self.file_compact_date1 = datetime.strptime(d_str, f)

                    d = file_date[4][:-3].split('-')[1]
                    f = "%Y-%m"
                    if int(d[4:6]) == 12:
                        d_str = str(int(d[0:4])+1) + '-' + str(1)
                    else:
                        d_str = str(d[0:4]) + '-' + str(int(d[4:6])+1)
                    self.file_compact_date2 = datetime.strptime(d_str, f) 


                    if str('_' + var + '_') in str(thing) and self.compact_date >= self.file_compact_date1\
                             and self.compact_date < self.file_compact_date2: 
 
                        data = Extract(str(self.three_month_path + thing), var, self.target_date, self.run_type)        
                        self.inst_list[var] = data
                            

        
        for thing in os.listdir(self.host_path):
            if 'constants' in thing:
                for var in meta_list:
                    data = Extract(str(self.host_path + thing), var, self.target_date, self.run_type)         
                    self.inst_list[var] = data



class main:  
        def __init__(self, input_date, model_run):
            self.model_run = model_run
            
            spec_date = input_date.split('-')[0] + '-' + input_date.split('-')[1] +\
                 '-' + input_date.split('-')[2] + 'T' + input_date.split('-')[-1] + ':00:00'

            Instance = Instance_Dict_Builder(input_date, self.model_run)
            input_date_instance = Instance_Functions(Instance.inst_list, input_date)
            
            try:

                # SET THESE FOR TRACKING PURPOSES
                self.up_hel = input_date_instance.up_hel[...]
                self.bunkers_uv = input_date_instance.bunkers_uv
                self.mean_wind = input_date_instance.mean_wind


            except AttributeError:
                pass                        

            # Grabs the reflectivity file directly. Not the most efficient considering it should be done with the input_date_instance.
	# But it works just fine.
            file_name = get_radar_file(input_date, self.model_run)

	# KEY LINE HERE: This calls the Pyartgridder to produce the PPyart object of reflectivity
	# Needed for tracking!
            self.grid1 = PyartGridder.wrf_make_grid(file_name,
            spec_date, self.model_run).pyart_grid

            self.longitude = PyartGridder.wrf_make_grid(file_name,
            spec_date, self.model_run).lonss 
            self.latitude = PyartGridder.wrf_make_grid(file_name,
            spec_date, self.model_run).latss 

            
            input_date = input_date.split('-')
            hour = float(input_date[-1])


# This class is called after the tracking is complete. 
# Updates UH, bunkers, mean wind, etc. in the tracking CSV file.
# Obviously these paths need to be relative.
class track_updates:
        
    def __init__(self, tracks_file, data_dict, start_date, sim, setup):
        self.setup = setup
        self.sim = sim   
        self.tracks_file = tracks_file

        self.data_dict = data_dict
        self.start_date = start_date

        self.get_uh()
        self.get_bunkers()
        self.future_locations()

    # Grabs the UH for the cell locations.
    def get_uh(self):

        self.tracks_file.to_csv('./tracks/nobunk_%s.csv' % (str(self.start_date)))
        read_file = pd.read_csv('./tracks/nobunk_%s.csv' % (str(self.start_date)))

        updraft_helicity = []

        for index, line in read_file.iterrows():
            try:
               date = str(line[2])
               date = str(date[:10]) + '-' + str(date[11:13])

               if int(date[11:13]) % 3 == 0:

                   lon, lat = float(line[5]), float(line[6])
                   i,j = get_closest_lat_lon(self.data_dict[date].latitude, 
                        self.data_dict[date].longitude, lat, lon)
                  
                   uh, ii,jj = get_nearby_max(self.data_dict[date].up_hel, i, j) 
                   updraft_helicity.append(self.data_dict[date].up_hel[ii,jj])

               else:
                   updraft_helicity.append('-')

            except ValueError:
                updraft_helicity.append('-')
         
        self.tracks_file['UH'] = updraft_helicity

    # Grabs the Bunkers motion vector for the cell locations. 
    def get_bunkers(self):

    read_file = pd.read_csv('./tracks/nobunk_%s.csv' % (str(self.start_date)))
    bunkers_u, bunkers_v, mean_u, mean_v = [], [], [], []

    for index, line in read_file.iterrows():
               
        date = str(line[2])
        hour = int(str(date[11:13]))
        date = str(date[:10]) + '-' + str(date[11:13])

        date = hour_check(line[2])
        
        date = str(date[:10]) + '-' + str(date[11:13])
        

        lon, lat = float(line[5]), float(line[6])

        i,j = get_closest_lat_lon(self.data_dict[date].latitude, 
                self.data_dict[date].longitude, lat, lon)


# Awkward section here. But this changes the points that are pulled for the Bunkers and Mean
# wind calculation if the cell is too close to the boundaries. 
        if i > 11: mindis_i = 10
        elif i < 11 and i > 1: mindis_i = i - 1
        else: mindis_i = 0

        if j > 11: mindis_j= 10
        elif j < 11 and j > 0: mindis_j = j - 1
        else: mindis_j = 0


        if 531 - j > 11: sm_j = 10
        elif 531 - j < 11 and 531 - j > 1: sm_j = 531 - j - 1
        else: sm_j = 0

        if 691 - i > 11: sm_i = 10
        elif 691 - i < 11 and 691 - i > 0: sm_i = 691 - i - 1
        else: sm_i = 0

        bunkers_u.append(np.ma.mean(np.ma.masked_invalid(self.data_dict[date].bunkers_uv[0][i-mindis_i:i+sm_i,j-mindis_j:j+sm_j])))

        bunkers_v.append(np.ma.mean(np.ma.masked_invalid(self.data_dict[date].bunkers_uv[1][i-mindis_i:i+sm_i,j-mindis_j:j+sm_j])))
        
        mean_u.append(np.ma.mean(np.ma.masked_invalid(self.data_dict[date].mean_wind[0][i-mindis_i:i+sm_i,j-mindis_j:j+sm_j])))

        mean_v.append(np.ma.mean(np.ma.masked_invalid(self.data_dict[date].mean_wind[1][i-mindis_i:i+sm_i,j-mindis_j:j+sm_j])))
                

        self.tracks_file['Mean_U'] = mean_u
        self.tracks_file['Mean_V'] = mean_v

        self.tracks_file['Bunkers_U'] = bunkers_u
        self.tracks_file['Bunkers_V'] = bunkers_v

    def future_locations(self):

            self.tracks_file.to_csv('./tracks/nopredict_%s.csv' % (self.start_date))
            
            read_file = open('./tracks/nopredict_%s.csv' % (self.start_date))
            predict_lat, predict_lon = [], []
            
            uids = {}              
            count = 0

            for l in read_file:
                line = l.split(',')

                if count == 0:
                    count += 1
                    continue        


                if line[0] == '0':
                    predicts = (line[6], line[5])
                    uids[str(line[1])] = {str(line[0]): [line[6],line[5],line[-2],line[-1]]}
                      
                else:

                    previous_hour = float(line[0]) - 1

                    try:
                        newline = {str(line[0]): [line[6],line[5],line[-2],line[-1]]}
                        uids[str(line[1])].update(newline)   

                        previous_bu = uids[str(line[1])][str(int(previous_hour))][2]
                        previous_bv = uids[str(line[1])][str(int(previous_hour))][3]
                        previous_lat = uids[str(line[1])][str(int(previous_hour))][0]
                        previous_lon = uids[str(line[1])][str(int(previous_hour))][1]
                        predicts = predict_location(float(previous_bu), float(previous_bv), 
                                                      float(previous_lat), float(previous_lon))
                    except KeyError:
                        uids[str(line[1])] = {str(line[0]): [line[6],line[5],line[-2],line[-1]]}
                        predicts = (line[6], line[5])

                    predict_x, predict_y = predicts[0], predicts[1]

                    predict_lat.append(predict_x)
                    predict_lon.append(predict_y)

            self.tracks_file['Predict_Lat'] = predict_lat
            self.tracks_file['Predict_Lon'] = predict_lon
            read_file.close()
            
            
            read_file = open('./tracks/nopredict_%s.csv' % (self.start_date))
            predict_lat, predict_lon = [], []
            
            uids = {}
            
            count = 0
            for l in read_file:
                 line = l.split(',')

                 if count == 0:
                        count += 1
                        continue        

                 if line[0] == '0':
                        predicts = (line[6], line[5])
                        uids[str(line[1])] = {str(line[0]): [line[6],line[5],line[-4],line[-3]]}
                        
                 else:
                 previous_hour = float(line[0]) - 1

                    try:
                        newline = {str(line[0]): [line[6],line[5],line[-4],line[-3]]}
                        uids[str(line[1])].update(newline)   

                        previous_mu = uids[str(line[1])][str(int(previous_hour))][2]
                        previous_mv = uids[str(line[1])][str(int(previous_hour))][3]
                        previous_lat = uids[str(line[1])][str(int(previous_hour))][0]
                        previous_lon = uids[str(line[1])][str(int(previous_hour))][1]

                        predicts = predict_location(float(previous_mu), float(previous_mv), 
                                                  float(previous_lat), float(previous_lon))
                    except KeyError:
                        uids[str(line[1])] = {str(line[0]): [line[6],line[5],line[-4],line[-3]]}
                        predicts = (line[6], line[5])

                    predict_x, predict_y = predicts[0], predicts[1]

                    predict_lat.append(predict_x)
                    predict_lon.append(predict_y)

                self.tracks_file['Mean_Predict_Lat'] = predict_lat
                self.tracks_file['Mean_Predict_Lon'] = predict_lon                


                self.tracks_file.to_csv('./tracks/rawtracks_%s.csv' % \
                  (str(self.start_date)[0:10]))





def get_params(setup_string):
   setup_split = setup_string.split('_')
   sm = int(setup_split[0][2:])
   fm = int(setup_split[1][2:])      
   dis = int(setup_split[2][2:])
   mfm = int(setup_split[3][2:])
   shift = int(setup_split[4][2:])
   return sm, fm, dis, mfm, shift

def looper(start, run_type, setup):
    full_date_dict = {}
    date_list = []
    
    input_date = start
    
    for i in range(0,25):

            date_list.append(input_date)
            grid1 = main(input_date, run_type)
            full_date_dict[input_date] = grid1


            year, month, day, hour = input_date.split('-')
            start = datetime(year=int(year),month=int(month),day=int(day),hour=int(hour))
            delta = timedelta(hours=1)
            offset = start + delta
            nexthour_time = offset.strftime('%Y-%m-%d-%H')
            input_date = nexthour_time
           


    grids_list = []                 
    for griddate in date_list:
            grids_list.append(full_date_dict[griddate].grid1)
    
    # Tracking Parameter Defaults
    FIELD_THRESH = 45
    ISO_THRESH = 8
    ISO_SMOOTH = 3
    MIN_SIZE = 20
    GS_ALT = 1500

    grids = (item for item in grids_list)
    
    SEARCH_MARGIN, FLOW_MARGIN, MAX_DISPARITY, MAX_FLOW_MAG, MAX_SHIFT_DISP = get_params(setup)

    tracks_obj = Cell_tracks()

    tracks_obj.params = {'FIELD_THRESH': FIELD_THRESH,
                   'MIN_SIZE': MIN_SIZE,
                   'SEARCH_MARGIN': SEARCH_MARGIN,
                   'FLOW_MARGIN': FLOW_MARGIN,
                   'MAX_FLOW_MAG': MAX_FLOW_MAG,
                   'MAX_DISPARITY': MAX_DISPARITY,
                   'MAX_SHIFT_DISP': MAX_SHIFT_DISP,
                   'ISO_THRESH': ISO_THRESH,
                   'ISO_SMOOTH': ISO_SMOOTH,
                   'GS_ALT': GS_ALT}


    tracks_obj.get_tracks(grids)
    tracks_file = tracks_obj.tracks        

    track_object = track_updates(tracks_file, full_date_dict, start, run_type, setup)




