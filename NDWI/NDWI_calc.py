#!/usr/bin/env python

"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Command line utility for calculating the NDWI for a input directory of bulk-Order
Landsat scenes from Earth Explorerself.

Need to run sensitivity tests on what is the best -pix_thres to use to most
accuratley display the distribution of water on the glaciers.

Since this script is most i/o bound (reading and writing rasters) I have not
attempted to parallelize any part of the processing step. Instead, I followed the
methods from - https://www.gis.usu.edu/~chrisg/python/2009/lectures/ospy_slides4.pdf
& http://geoexamples.blogspot.com/2013/06/gdal-performance-raster-classification.html
The raster is read in and written by a user specified block size (default 256 x 256)
instead of the deafualt Landsat blocksize of 9356 x 1. This is suppoded to be the
far more efficent way to read, process and write the data the reading the entire
raster in one sitting, processing, and writing in one sitting..

Another option may be to parrallelize the decompression of the .tar.gz Landsat
scenes. This will present a significant slow down if it is only computed by one
core. This will especially be the case when the script is executed on an direcotry
containg 50-100 Landsat scenes to be processed.

TODO:
	- Attemot to .tar.gz files in parallel, using pigz?
		- need to find out if really better to do in parallel, will depend on how
		  much of the process is cpu bound and how much is i/o bound.
	- Do I need to do Top of Atmosphere Reflecatance (TOA) corrections for the
	  images since they are level-1 or should I just download higher level scenes
	  (ie. Level-2 or ready for processing)
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

import os
import osr
import sys
import gdal
import argparse
import subprocess
import numpy as np
from argparse import RawTextHelpFormatter

def espg_code(gdal_ds):

    proj = osr.SpatialReference(wkt=gdal_ds.GetProjection())
    espg = int(proj.GetAttrValue('AUTHORITY',1))

    return espg

def calc_ndwi(Green, NIR, QA, NDWI_empty):
    valid_pixels = (QA !=1) & (Green + NIR != 0)

    NDWI_empty[:] = -9999
    Green = Green[valid_pixels].astype(np.float64)
    NIR = NIR[valid_pixels].astype(np.float64)
    NDWI_empty[valid_pixels] = (Green - NIR) / (Green + NIR)

    return NDWI_empty

def ndwi_mask(NDWI,pix_thresh=0.15):
    valid_pixels = (NDWI != -9999) & (NDWI >= pix_thresh)

    NDWI_mask= np.empty_like(NDWI, np.float64)
    NDWI_mask[:] = 0

    NDWI = NDWI[valid_pixels].astype(np.float64)

    NDWI_mask[valid_pixels] = 1

    return NDWI_mask


def reproject(raster_fn,t_srs=None):
    '''
        -   Right have it defaulting to Polar Stereo and manually inputting the parameters. Could and should change to use the EPSG code to alloow for
            other projections.

        -   Ideally wouldn't need to do this on the command line but due to errors with GDAL_DATA variable this is the sutiable alternative for right now
    '''
    no_f_type = os.path.splitext(raster_fn)[0]
    PS_fn = no_f_type + 'PS.TIF'

    if os.path.exists(PS_fn) == False:
        print ('Raster will be reprojected with Polar Stereo (ESPG: 3413)')
        #warp_cmd = "gdalwarp -t_srs '+proj=stere +lat_ts=70 +lat_0=90 +lon_0=-45+y_0=0 +x_0=0 +k=1 +datum=WGS84 +units=m'" + ' ' + raster_fn + ' ' + PS_fn
        warp_cmd = 'gdalwarp -t_srs "EPSG:' + str(t_srs) + '"'+ ' ' + raster_fn + ' ' + PS_fn
        subprocess.call(warp_cmd, shell=True)
    return PS_fn

def read_n_write(green_fn, nir_fn, qa_fn, t_srs=None, x_block_size = 256, y_block_size = 256):


    green_ds = gdal.Open(green_fn)
    nir_ds = gdal.Open(nir_fn)
    qa_ds = gdal.Open(qa_fn)

    #All data sets should have the same projection so only test the projection of the first one

    if t_srs == None:
        green_band = green_ds.GetRasterBand(1)
        nir_band = nir_ds.GetRasterBand(1)
        qa_band = qa_ds.GetRasterBand(1)
    else:
        espg = espg_code(green_ds)
        if t_srs != espg:
            green_rpjct = reproject(green_fn,t_srs)
            nir_rpject = reproject(nir_fn, t_srs)
            qa_rpject = reproject(qa_fn, t_srs)

            green_ds = gdal.Open(green_rpjct)
            nir_ds = gdal.Open(nir_rpject)
            qa_ds = gdal.Open(qa_rpject)

            green_band = green_ds.GetRasterBand(1)
            #green_band.astype(np.float64)
            nir_band = nir_ds.GetRasterBand(1)
            #nir_band.astype(np.float64)
            qa_band = qa_ds.GetRasterBand(1)
            #nir_band.astype(np.float64)

    # All input raster are of the same size and the output raster will be of the same size aswell
    # so only getting the dimensions of the raster once to be effeicent
    rows = green_ds.RasterYSize
    cols = green_ds.RasterXSize

    xbs = x_block_size
    ybs = y_block_size

    # only here temporalily to make sure the array is 0 and 1

    whole_array = np.empty([rows,cols])

    #sets up the output raster that will be written to block by block
    #Determine filename, again only using the first input since info should be the same for all the rasters
    parent_dir, fn = os.path.split(green_fn)
    mask_fn = parent_dir + '/'+'_'.join(fn.split('_')[:-1])+'_NDWI_mask.TIF'
    print(mask_fn)
    #create the data set, define transform and projection
    outdriver = gdal.GetDriverByName("GTiff")
    NDWI_ds = outdriver.Create(mask_fn, cols, rows, 1, gdal.GDT_Int16)
    NDWI_ds.SetGeoTransform(green_ds.GetGeoTransform())
    NDWI_ds.SetProjection(green_ds.GetProjection())

    # finding the index of the subset block of the raster to be read in
    for i in range(0, rows, ybs):
        if i + ybs < rows:
            numRows = ybs
        else:
            numRows = rows - i
        for j in range(0, cols, xbs):
            if j + xbs < cols:
                numCols = xbs
            else:
                numCols = cols - j

            green_array = green_ds.ReadAsArray(j, i, numCols, numRows)
            nir_array = nir_ds.ReadAsArray(j, i, numCols, numRows)
            qa_array = qa_ds.ReadAsArray(j, i, numCols, numRows)
            NDWI_array = np.empty([numRows,numCols], np.float64)

            NDWI = calc_ndwi(green_array, nir_array, qa_array, NDWI_array)
            NDWI_mask = ndwi_mask(NDWI, pix_thresh=0.15)
            whole_array[i:(i + numRows), j:(j + numCols)] = NDWI_mask
            #This is where the script is getting hung up

            NDWI_ds.GetRasterBand(1).SetNoDataValue(-9999)
            NDWI_ds.GetRasterBand(1).WriteArray(NDWI_mask,j,i)

    NDWI_ds= None

    return whole_array


def getparser():
    description = ('Command line utility to calculate the Normailzed Differenced Water Index (NDWI) for a bulk order directory from earth Explorer. ')
    parser = argparse.ArgumentParser(description=description,formatter_class=RawTextHelpFormatter)
    parser.add_argument('src_dir',type=str,help='The source directory. Should be all of the available Landsat Scenes for a glacier. Should be a Bulk-Order directory from Earth Explorer (Ex. "Bulk_Order_930842")')
    parser.add_argument('-t_srs',type=str, default = None, help='The target projection is the Landsat Scenes need to be reprojected')
    parser.add_argument('-pix_thres',type=float, default = 0.15, help = 'minimum value for the NDWI pixel value to be considered valid water.')
    parser.add_argument('-m',action='store_false',help='Whether the output should be a mask. Add flag to command line input to output mask.')
    parser.add_argument('-x_block_size', default = None, type = int, help = 'block size in the x direction')
    parser.add_argument('-y_block_size', default = None, type = int, help = 'block size in the y direction')
    return parser

def main():
    parser = getparser()
    args = parser.parse_args()

    src_dir = args.src_dir
    t_srs = args.t_srs
    pix_thres = args.pix_thres

    #test that if either x_block_size or y_block_size is defined the other is also defined
    if args.x_block_size != None and args.y_block_size == None:
            parser.error('\n\n-x_block_size also requires -y_block_size \n\n')
    elif args.y_block_size != None and args.x_block_size == None:
        parser.error('\n\n -y_block_size also requires -x_block_size \n\n')
    elif args.y_block_size == None and args.x_block_size == None:
        x_block_size = 256
        y_block_size = 256
    else:
        x_block_size = args.x_block_size
        y_block_size = args.y_block_size




    cwd = os.getcwd()

    for root,dirs,files in os.walk(src_dir):
        '''
        Needs to be sorted out. The inital walk through the directory with the files still zipped show up as file and not as dirs. But once they
        unzipped they will be show up as dirs which will complicate the os.walk process since they will be located in seperate list.
        '''

        # This is the inital walk for when the files are still zipped.
        # Could be a  good idea to parallelize this. Multithreading this and giving each core a direcotry to unzip might not be the worst idea
        # especially when you will have 100+ scenes per direcotry of each glacier.
        for line, fn in enumerate(files):
            #compile a list of all of the scene directories and their file paths to itterate through
            if fn.split('_')[0][2:] == '08':
                if fn.endswith('.tar.gz') == True and os.path.isdir(root + '/'+ fn.split('.')[0]) == False:
					unzipped_dir = root + '/' + fn.split('.')[0]
                    os.mkdir(unzipped_dir)
                    unzip_cmd = 'tar -xvzf ' + root+ '/' + fn + ' -C ' + unzipped_dir
                    subprocess.call(unzip_cmd,shell=True)

					#more effiecent way to unzipp and the imdeately do the necessary calculations. Otherwise would have created bottleneck
					# as the 50-100 direcotries unzipp. 
					for i,band in enumerate(os.listdir(unzipped_dir)):
						if band.endswith('B3.TIF') == True:
							green_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                        elif band.endswith('B5.TIF') == True:
                            nir_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                        elif band.endswith('BQA.TIF') == True:
                            qa_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]

                        whole_array = read_n_write(green_fn, nir_fn, qa_fn, t_srs=t_srs, x_block_size = x_block_size, y_block_size = x_block_size)


                else:
                    # This is the second walk for once the files are unzipped and appear as dirs.
                    for line, dn in enumerate(dirs):
                        #Might need to add more conditionals here based on what the naming convention is for all of Landsat 8
                        if dn.startswith('LC') == True:
                            bands = os.listdir(root + '/' +dn)
                            for i,band in enumerate(bands):
                                if band.endswith('B3.TIF') == True:
                                    green_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                                elif band.endswith('B5.TIF') == True:
                                    nir_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]
                                elif band.endswith('BQA.TIF') == True:
                                    qa_fn = cwd + '/' + root + '/' + dn + '/' + bands[i]


                            whole_array = read_n_write(green_fn, nir_fn, qa_fn, t_srs=t_srs, x_block_size = x_block_size, y_block_size = x_block_size)

if __name__ == '__main__':
    main()
