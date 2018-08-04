#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  2 09:43:03 2018

@author: andrewnolan
"""
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
import geopandas as gpd
from shapely.geometry import Point, Polygon

description = ('Write .txt file of PGC DEM IDs of DEMs to request. Data is pulled based on central coordinates of terminus. \n \ndependecies: \n \t-shaply \n \t-geopandas')

def getparser():
    parser = argparse.ArgumentParser(description=description,formatter_class=RawTextHelpFormatter)
    parser.add_argument('shp_fn',type=str,help='Input Shp files name and filepath')
    parser.add_argument('term_coords',type=str, help='terminus coordinates (Deimcal Degrees) for the glacier of interest -format:(lat,lon)')
    parser.add_argument('-out_fn',type=str,help='filepath and name for the .txt file to be written')
    return parser
def main():
    parser = getparser()
    args = parser.parse_args()

    #Check if shp file exists
    shp_fn = args.shp_fn
    if not os.path.exists(shp_fn):
        sys.exit("Unable to find shp_fn: %s" % shp_fn)

    #Read shp file into geopandas df
    data = gpd.read_file(shp_fn)

    #take terminus position and create a shply point geometry object
    term_coords = args.term_coords
    lat = float(term_coords.split(',')[0])
    lon = float(term_coords.split(',')[1])
    term_point = Point(lon,lat)

    out_fn = args.out_fn
    if out_fn is None:
        out_fn = os.getcwd()+'/file_names.txt'


    #OLD METHOD
    #create polygon from centroid point
    #p = Point(data['CENT_LAT'][i],data['CENT_LON'][i]).buffer(data['Shape_Area'][i])

    #find the filepaths that the coordinate lies within
    PAIRNAME_list = []
    for i in range(len(data)):
        #itterate through each polygon to test if the point is within it and test wether DEM is along track
        if term_point.within(data['geometry'][i]) and data['DEM_ID'][i][0:2]=='WV':
            dem_PAIRNAME = str(data['PAIRNAME'][i])
            PAIRNAME_list.append(dem_PAIRNAME)


    PAIRNAME_set = set(PAIRNAME_list)
    output_file =open(out_fn,'w')
    for i in PAIRNAME_set:
        output_file.write('{0:}''\n'.format(i))
    output_file.close()

if __name__ == '__main__':
    main()
