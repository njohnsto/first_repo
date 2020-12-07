"""
This attached script is provided as a tool and NOT as a RaySearch endorsed script for
clinical use.  The use of this script in its whole or in part is done so
without any guarantees of accuracy or expected outcomes. Please verify all results. Further,
per the Raystation Instructions for Use, ALL scripts MUST be verified by the user prior to
clinical use.

Script name : CreateROImask_wMARGINfunc_CalculateHUforROI_2.py
Author      : Tibor Marusic, RaySearch Laboratories
Description : Script finds the HU values for a specific ROI and returns the average
			  and standard deviation values for it.
			 
NOTE		: - It currently only works for ROIs of type Contour , not MESH
			  - It assumes that all Slice distances are constant
			  - NOT TESTED FOR FFP or any other otrientation except FFS
			  - Does not yield the exact values that RS does, probably because
				the algorithm for choosing voxels (or even interpolating them for 
				partially covered voxels) is different in RS compared to the script

Prerequisite: Needs to have Numpy and OpenCV modules installed and is only
			  tested for Python 3.6
"""

# import matplotlib.pyplot as plt #Only for plotting purposes
import numpy as np
import cv2
from connect import *


'''
DEFINING FUNCTIONS
'''

def convert_ROI_to_voxel_FOR(ROI_coord, corner = np.array([0,0,0]), voxel_dim = np.array([1,1,1]) ):
	'''
	Function to convert ROI coordniates into voxel FOR, meaning that coord (0,0,0) is the start
	Please note that voxel_dim = (PixelSize_x, PixelSize_y, Slice_size_z)
	and corner = (corner_x, corner_y, corner_z)
	
	NOTE:
	THIS ASSUMES THAT THE CT SLICE distance IS CONSTANT !!!!
	'''
	coord = (ROI_coord - corner) / voxel_dim 
	
	return coord


def vdx_slice_data(data, x=None, y=None, z=None):
	'''
	returns all the data in the same slice where x, y or z= a value
	does not work if x, y and z have values at the same time
	have not yet added restriction for choosing only 1 dimention , could do 0=x, 1=y and 2=z , maby easiest solution.
	'''

	if x != None:
		list_of_points=[i for i in range(len(data)) if data[i,0]==x]
		return data[list_of_points]
	
	if y != None:
		list_of_points=[i for i in range(len(data)) if data[i,1]==y]
		return data[list_of_points]
	
	if z != None:
		list_of_points=[i for i in range(len(data)) if data[i,2]==z]
		return data[list_of_points]

def vdx_volume_mask(voi_data, dimx, dimy, dimz, margin=None):
	'''
	Attempting to generate a 3d mask of the ROI (voi).
	voi_data: Are the ROI coordinates given as an numpy array of (x,y,z) coordinates
	dimx, dimy and dimz: are the dimetntions of the mask.
	margin: You can add a distance margin to the mask to include more voxels surrounding the ROI
			If left empty there will be no margin made
	'''
	mask_size_voi=np.array([dimx, dimy, dimz])#can be changed to any number, example 120 (=120x120x120 pixels)
	
	mask_voi=np.zeros((dimx,dimy,dimz))
	
	voi_data[:,2] = np.round(voi_data[:,2])
	
	
	#Generating mask slices, also some overlap algebra is used to try to exclude the effects of having ROI coordinates differ from CT resolution
	for z_nr in range(dimz):
		polygons=vdx_slice_data(voi_data, z=z_nr)
		if len(polygons) != 0:
			mask_temp_nn=cv2.fillPoly(np.array(mask_voi[:,:,z_nr]), pts = [np.array(np.round(polygons[:,0:2]+np.array([-0.5,-0.5]))).astype(int)], color = (1,1,1))
			
			mask_temp_np=cv2.fillPoly(np.array(mask_voi[:,:,z_nr]), pts = [np.array(np.round(polygons[:,0:2]+np.array([-0.5,0.5]))).astype(int)], color = (1,1,1))
			
			mask_temp_pn=cv2.fillPoly(np.array(mask_voi[:,:,z_nr]), pts = [np.array(np.round(polygons[:,0:2]+np.array([0.5,-0.5]))).astype(int)], color = (1,1,1))
			
			mask_temp_pp=cv2.fillPoly(np.array(mask_voi[:,:,z_nr]), pts = [np.array(np.round(polygons[:,0:2]+np.array([0.5,0.5]))).astype(int)], color = (1,1,1))
			
			
			mask_temp=mask_temp_nn*mask_temp_np*mask_temp_pn*mask_temp_pp
			mask_voi[:,:,z_nr]=np.flip(np.rot90(mask_temp), axis=0)
	
	
	#ADDING MARGIN
	if margin != None:
		voi_coord = np.copy(np.transpose(np.nonzero(mask_voi)))
		for i in voi_coord:
			for x in np.arange(i[0]-margin, i[0]+margin+1):
				for y in np.arange(i[1]-margin, i[1]+margin+1):
					for z in np.arange(i[2]-margin, i[2]+margin+1):
						mask_voi[x, y, z] = 1
						
	
	# cv2.imshow(" ", newimage) # this is only to show the image
	# cv2.waitKey()
			

	return mask_voi.astype(bool)
	
	
'''
CODE STARTS HERE
'''
#######################
### Choose ROI_NAME ###
#######################

ROI_NAME = 'PTV'

### Retrieving data from the current examination ###
exam = get_current('Examination')
case = get_current('Case')

dimx = exam.Series[0].ImageStack.NrPixels.x
dimy = exam.Series[0].ImageStack.NrPixels.y
dimz = len(exam.Series[0].ImageStack.SlicePositions)

corner = exam.Series[0].ImageStack.Corner
corner = np.array([corner.x, corner.y, corner.z])

pixel_size = exam.Series[0].ImageStack.PixelSize
slice_postion_array = exam.Series[0].ImageStack.SlicePositions
voxel_dim = np.array([pixel_size.x,pixel_size.y,slice_postion_array[1]-slice_postion_array[0]])

rescale_intercept = exam.Series[0].ImageStack.ConversionParameters.RescaleIntercept
rescale_slope = exam.Series[0].ImageStack.ConversionParameters.RescaleSlope
pixel_representation = exam.Series[0].ImageStack.ConversionParameters.PixelRepresentation


# ### Converting pixel data from a byte array to 3D intensity cube ###
# ### NOTE: THIS PART IS WHAT MAKES THE CODE RUN SLOW !!!!         ###
# pixel_data = exam.Series[0].ImageStack.PixelData
# array = []
# for i in range(int(len(pixel_data)/2)):
	# array.append(pixel_data[i*2+1]*256+pixel_data[i*2])
# HU_data = np.reshape(array,(dimx,dimy,dimz), order = 'F')   #np.rot90(, 3) -- No need apperently


'''
NEW WAY - MUCH FASTER
Also with support for signed CT values
'''
pixel_data = exam.Series[0].ImageStack.PixelData
pixel_data = pixel_data.astype(np.float)
length = len(pixel_data)
evens = np.arange(0, length, 2, dtype=np.int)
odds = np.arange(1, length, 2, dtype=np.int)
if pixel_representation == 0:
	array = pixel_data[evens]+pixel_data[odds]*256
else:
	array = pixel_data[evens]+pixel_data[odds]*256
	array = array.astype(np.int16)
HU_data = np.reshape(array,(dimx,dimy,dimz), order = 'F')

### Converting ROI coordinates to (x,y,z)- array ###
ROI_coord_obj = case.PatientModel.StructureSets[exam.Name].RoiGeometries[ROI_NAME].PrimaryShape.Contours
ROI_coord = []
for i in ROI_coord_obj:
	for j in i:
		ROI_coord.append([j.x,j.y,j.z])
ROI_coord = np.asarray(ROI_coord)


### Converting ROI coordinates to CT frame of reference ###
ROI_coord_voxel_FOR = convert_ROI_to_voxel_FOR(ROI_coord, corner, voxel_dim)


### Creating a mask that corresponds to only ROI voxels ###
ROI_mask = vdx_volume_mask(ROI_coord_voxel_FOR, dimx, dimy, dimz, margin=None)


### Calculating average and STD whilst converting to HU units ###
average = np.average(HU_data[ROI_mask])*rescale_slope+rescale_intercept 
std = np.std(HU_data[ROI_mask]*rescale_slope+rescale_intercept)

print('Average:', average)
print('\nStd:', std)


'''
Plotting purposes
'''
# cv2.imshow(" ", HU_data[:,:,10]/256) 
# cv2.waitKey()
# plt.imshow(HU_data[:,:,10],cmap='gray',interpolation='none')
# plt.show()



