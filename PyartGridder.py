from cartopy.feature import NaturalEarthFeature
import cartopy.crs as crs
import netCDF4 as nc4
from pyart import core, io
import numpy as np
import matplotlib.pyplot as plt
from tint import Cell_tracks, animate
from pprint import pprint
import re 

def get_point(latpoint, lonpoint, lats, lons):
	a = abs( lats-latpoint ) + abs( lons-lonpoint )
	i,j = np.unravel_index(a.argmin(), a.shape)
	return i,j

# Takes a 2D grid of reflectivity and makes a Py-Art grid object that can be input into
# TiNT readily.

# Some of this is can be 
class wrf_make_grid:

	def __init__(self, filename, date, run_type):

		reformdate = [str(int(item)) for item in re.split('-+|T+|:', str(date))]
		self.date = reformdate[0] + '-' + reformdate[1] + '-' + reformdate[2] + 'T' + reformdate[3] + ':' + '00' + ':' + '00'
		year = reformdate[0]

		self.run_type = run_type
		self.data = nc4.Dataset('%s' % filename, 'r')
		
        # 2D lat-lon projection. 
		self.latss = self.data.variables['XLAT']
		self.lonss = self.data.variables['XLONG']

        middle_x, middle_y = len(self.latss[:,0]) / 2., len(self.latss[0,:]) / 2.

		self.x, self.y = core.geographic_to_cartesian(self.lonss, self.latss, {'lat_0':self.latss[middle_x, middle_y],'lon_0':self.lonss[middle_x, middle_y], 'proj':'lcc'})


		self.xpoints, self.ypoints = self.x[middle_x,:], self.y[:,middle_y]
	
		self.zpoints = np.array([0])

		self.radar_longitude = {'data': [self.lonss[middle_x, middle_y]]}
		self.radar_latitude = {'data': [self.latss[middle_x, middle_y]]}
		self.radar_altitude = {'data': [1500]}

		
		self.proj = {'lat_0':self.latss[middle_x, middle_y],\  
            'lon_0':self.lonss[middle_x, middle_y], 'proj': 'lcc'}
		self.origin_latitude = {'standard_name':'latitude', 'data':self.latss[middle_x, middle_y]}
		self.origin_longitude = {'standard_name':'longitude', \
            'data':self.lonss[middle_x, middle_y]}
		self.x = {'data': self.xpoints}
		self.y = {'data': self.ypoints}
		self.z = {'data': self.zpoints}	
	
		self.get_time_index()


		self.time1 = {'standard_name': 'time', 'calendar': 'gregorian',
		'long_name': 'Time of grid', 'units':'seconds since time %sZ' % 
		self.date, 'data': np.array([self.date])}
		
		self.fields = {}

		datas = np.ma.masked_where(self.data.
		variables['REFL_10CM'][self.time_index,...] < 10, 
		self.data.variables['REFL_10CM'][self.time_index,...])

		datas_ = datas[np.newaxis,...]

		self.field1 = {'reflectivity': {'_FillValue': -9999.0, 
		'data': datas_}}

		self.origin_altitude = {'standard_name':'altitude', 'data':[0.]}

		self.metadata = {'Conventions': 'Simulated WRF 10cm'}


		self.grid = core.Grid(self.time1, self.field1, self.metadata, self.origin_latitude, self.origin_longitude, self.origin_altitude, self.x,self.y,self.z,projection=self.proj,radar_latitude=self.radar_latitude, radar_longitude=self.radar_longitude, radar_altitude=self.radar_altitude)

		self.grid.point_longitude['data'] = np.array(self.lonss)
		self.grid.point_latitude['data'] = np.array(self.latss)

		self.get_cell_locations(self.grid)


	def get_cell_locations(self,grid1):	
		self.pyart_grid = grid1
		grids1 = [self.pyart_grid]

	def modeltime_to_datetime(self):
		start = datetime(year=1901,month=1,day=1,hour=0)
		delta = timedelta(hours=1) * int(self.modeltimes)
		offset = start + delta
		single_time = offset.strftime('%Y-%m-%dT%H:%M:%S')
		return single_time


	def get_time_index(self):
		for i, key in enumerate(self.data.variables['Times']):
			strdate = ''
			for item in key:
				item = item.decode('UTF-8')
				strdate = strdate + item
			strdate = strdate.replace('_','T')
			strdate = [str(int(item)) for item in re.split('-+|T+|:', str(strdate))]
			strdate = strdate[0] + '-' + strdate[1] + '-' + strdate[2] + 'T' + strdate[3] + ':' + '00' + ':' + '00'


			if str(self.date) in str(strdate): 
				self.time_index = i


