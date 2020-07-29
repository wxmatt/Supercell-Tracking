import multiprocessing as mp
import re 
from datetime import datetime, timedelta
import csv
import matplotlib
import pandas as pd
import sys
import os
import math
import main_class
from math import sin, cos, sqrt, atan2
from tools import hour_check, predict_location
import matplotlib.pyplot as plt
from cartopy.feature import NaturalEarthFeature
import cartopy.crs as crs
import glob

# Function  to predict new lat-lon location based on previous location and velocity. 
def predict_location(bu, bv, start_lat, start_lon, direction):
   # Position, decimal degrees
   lat = start_lat
   lon = start_lon

   # Earths radius, sphere
   R = 6378137.

   # Offsets in meters
   if 'forward' in direction:
      dlat = bv * 60. * 60. 
      dlon = bu * 60. * 60.

   elif 'backward' in direction:
      dlat = -1. * bv * 60. * 60. 
      dlon = -1. * bu * 60. * 60.

   # Coordinate offsets in radians
   dLat = dlat / R
   dLon = dlon / (R * np.cos(np.pi * start_lat / 180.))

   # Offset Position, decimal degrees
   end_lat = lat + dLat * 180. / np.pi
   end_lon = lon + dLon * 180. / np.pi 

   return [end_lat, end_lon]
	
def distance(lat1, lat2, lon1, lon2):
   R = 6373.0

   dlon = math.radians(lon2 - lon1)
   dlat = math.radians(lat2 - lat1)

   a = (sin(dlat/2.))**2 + cos(math.radians(lat1)) * cos(math.radians(lat2)) * (sin(dlon/2.))**2.

   c = 2. * atan2(sqrt(a), sqrt(1-a))
   dist = R * c

   return dist   


# Produces the dates for which to combine tracking data over. 
# This currently runs from March - June 2001 to 2013
def make_date_list():

   start_dates = []

   years = range(2001,2014)   
   for count, year in enumerate(years):

      date = datetime(year, 3, 1, 0, 0,0)
      str_date = str(date)
      str_date = re.split('-| |:', str_date)
      str_date = str_date[0]+'-'+str_date[1]+'-'+str_date[2]

      # Loops from over days from date (122 right now for March - June)
      for i in range(122): 
          str_date = str(date)
          str_date = re.split('-| |:', str_date)
          str_date = str_date[0]+'-'+str_date[1]+'-'+str_date[2]
          start_dates.append(str_date)
          date += timedelta(days=1)

   return start_dates

# Updating the UIDs in the tracks_with_UH tracking files each day. 
def update_uids(run_type, dirr):

   file_start = './tracks/rawtracks_'

   date_list = make_date_list()

   for i, filedate in enumerate(date_list):
      file_name = file_start + date_list[i] + '.csv'

      try:
         
         current_df = pd.read_csv(file_name)
      
         # Add a large value to the UID to avoid duplications. This way each cell has unique UID throughout the time period. 
         try:
            current_df['uid'] = current_df['uid'] + ((1000000 + i * 10000) * i)
         except KeyError:
            pass
      
         current_df.to_csv('./tracks/tracks_with_added_uid_%s.csv' % (date_list[i]))
      
      except (IOError, pd.errors.EmptyDataError) as err:
         print(err) 

   
# Connect the cells that cross over from each daily tracking file. Unlikely since the tracking files
# go from 12z to 12z
def cross_overs(run_type, dirr):

   file_start = './tracks/tracks_with_added_uid_'

   date_list = make_date_list()

   for i, filedate in enumerate(date_list):
      try:
         file_name = file_start + date_list[i+1] + '.csv'
         
      except IndexError:
         break

      try:
         current_df = pd.read_csv(file_name)

         # Add value to the UID to avoid duplications. This way each cell has unique UID throughout the time period. 

         previous_df = pd.read_csv(file_start + date_list[i] + '.csv')
         try:
            try:
               if np.array(previous_df.loc[(previous_df['scan'] == 24),'lat']) == np.array(current_df.loc[(current_df['scan'] == 0), 'lat']):
                  vals = np.array(previous_df.loc[(previous_df['scan'] == 24), 'uid'])
                  current_df.loc[(current_df['scan'] == 0), 'uid'] = vals

            except ValueError:
               pass

            try:
                new_concat = pd.concat([new_concat,current_df]).drop_duplicates(['time','lat','lon'],keep='first').sort_values('time')
      
            except UnboundLocalError:
               new_concat = pd.concat([previous_df,current_df]).drop_duplicates(['time','lat','lon'],keep='first').sort_values('time')

         except KeyError:
            pass
      except IOError:
         pass   
   
   
   headers = list(new_concat.dtypes.index)

   # This is the final output file with the entire tracking dataset. 
   new_concat.to_csv('./tracks/tracks_concat.csv', columns=headers[1:])
   


import sys
run_type_dict = {}

# The directory with all the tracking stuff. 
directories = sorted(glob.glob('./tracks'))


# Opens the full concatenated tracks files, creates a dictionary of each supercell.
# Checks each supercell dictionary to confirm if the cell is a supercell
# Writes, total_tracks.csv as the final csv with all confirmed supercells

for full_directory in directories:

   directory = full_directory.split('/')[-1].split('-')[-1]
   
   if sys.argv[1] == '1':
      update_uids(directory)
      cross_overs(directory)
   
   # -------------------------------------------------------------------------------------------- #
   # -------------------------------------------------------------------------------------------- #
   # -------------------------------------------------------------------------------------------- #
   
   def check_uid(uid, dictionary):
      if uid in dictionary.keys():
         return True
      else:
         return False

   # -------------------------------------------------------------------------------------------- 

   cells = {}
   
   read_file = open('./tracks/tracks_concat.csv' % (directory))

   for index, line in enumerate(read_file):
      things = line.split(',')[1:]

      if index == 0:
         continue

      try:
         cells[str(things[1])]['locations'].extend([float(things[6]),float(things[5])])
         cells[str(things[1])]['type'].extend([things[11]])
         cells[str(things[1])]['reflect'].extend([float(things[9])])
         cells[str(things[1])]['area'].extend([float(things[7])])
         #print(things[12])

         if things[12] == '-': 
            uh = -999
         else:
            uh = float(things[12])

         cells[str(things[1])]['uh'].extend([uh])
         cells[str(things[1])]['mean_u'].extend([float(things[13])])
         cells[str(things[1])]['mean_v'].extend([float(things[14])])
         cells[str(things[1])]['bunkers_u'].extend([float(things[15])])
         cells[str(things[1])]['bunkers_v'].extend([float(things[16])])
         cells[str(things[1])]['time'].extend([str(things[2])])
      


      except KeyError:
         if things[12] == '-': 
            uh = -999
         else:
            uh = float(things[12])

         cells[str(things[1])] = {'locations': [float(things[6]), float(things[5])],\
            'type': [things[11]],'area': [float(things[7])], 'reflect': [float(things[9])], 'uh': [uh], \
            'mean_u':[float(things[13])], 'mean_v':[float(things[14])], 'bunkers_u':[float(things[15])],\
            'bunkers_v':[float(things[16])], 'time': [str(things[2])]}



   # Get the bunkers and mean wind predicted locations.
   for uid in cells.keys():
      uid = str(uid)
      
      if len(cells[uid]['time']) == 1:
         cells[uid]['bunkers_locations'] = cells[uid]['locations'][:2]
         cells[uid]['mean_locations'] = cells[uid]['locations'][:2]

      # def predict_location(bu, bv, start_lat, start_lon, direction):

      elif len(cells[uid]['time']) > 1:
         b_u = cells[uid]['bunkers_u']
         b_v = cells[uid]['bunkers_v']
         m_u = cells[uid]['mean_u']
         m_v = cells[uid]['mean_v']

         lats, lons = cells[uid]['locations'][0::2], cells[uid]['locations'][1::2]

         for loc in range(len(lats)):

            if loc == 0:
               cells[uid]['bunkers_locations'] = cells[uid]['locations'][:2]
               cells[uid]['mean_locations'] = cells[uid]['locations'][:2]

            else:                 
               hour = int(datetime.strptime(str(cells[uid]['time'][0]), '%Y-%m-%d %H:%M:%S').hour)
               
               b_u_mean, b_v_mean = np.mean(b_u[loc-1:loc+1]), np.mean(b_v[loc-1:loc+1])
               m_u_mean, m_v_mean = np.mean(m_u[loc-1:loc+1]), np.mean(m_v[loc-1:loc+1])

               cells[uid]['bunkers_locations'].extend(predict_location(b_u_mean, b_v_mean, lats[loc-1], lons[loc-1], 'forward'))
               cells[uid]['mean_locations'].extend(predict_location(m_u_mean, m_v_mean, lats[loc-1], lons[loc-1], 'forward'))


   # Now checking if supercell!
   uid_list, uhs, all_, all_bunks, bunks_highuh, all_uh = [], [], [], [], [], []
   i = 1
   for uid in cells.keys():
   
      uid = str(uid)

      bunkers_dist, mw_dist = [], []

      for q, num in enumerate(np.array(cells[uid]['locations'][0::2])):
         bunkers_dist.append(distance(cells[uid]['locations'][0::2][q], cells[uid]['bunkers_locations'][0::2][q], cells[uid]['locations'][1::2][q], cells[uid]['bunkers_locations'][1::2][q]))
        

      for q, num in enumerate(np.array(cells[uid]['locations'][0::2])):
         mw_dist.append(distance(cells[uid]['locations'][0::2][q], cells[uid]['mean_locations'][0::2][q], cells[uid]['locations'][1::2][q], cells[uid]['mean_locations'][1::2][q]))     
     
   
      mean_bunk, mean_mw = np.mean(bunkers_dist[1:]), np.mean(mw_dist[1:])
      

      if any(item > 45 for item in cells[str(uid)]['reflect']) and any(item == 'True' for item in cells[str(uid)]['type']):

         if (any(item >= 75 for item in cells[str(uid)]['uh'])) or ((mean_bunk < mean_mw) and mean_bunk < 30 and (any(item < 75 and item >= 25 for item in cells[str(uid)]['uh'])) and len(cells[str(uid)]['uh']) > 1 and (all(item < 75 for item in cells[str(uid)]['uh']))):

            uid_list.append(uid)  

   # Save the uid list for the current supercell configuration. Can change this obviously. 
   np.save('./UIDs/final/%s_uids_%s' % (run_type, 'All'), np.array(uid_list))
      
   with open('./tracks/total_tracks.csv', 'w') as csv_file:
      writer = csv.writer(csv_file)

      # CHANGE THIS LINE TO USE WHATEVER UID LIST YOU WANT.
      for key in uid_list:
            writer.writerow([key, cells[key]])

















   

