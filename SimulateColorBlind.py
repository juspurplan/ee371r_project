#!/usr/bin/python

# std python imports
import os.path

# other imports
import numpy
import scipy

import matplotlib
import matplotlib.pyplot as pyplot

import skimage
import skimage.color
import skimage.io

import click


# source for empirical studies on finding copunctal points. (Intersection points of 
# confusion lines for dichromatic viewers
# http://www.sciencedirect.com/science/article/pii/0042698996000892
R_COPUNCTAL = (.749, .251)
R_BLIND_ANGLE = scipy.pi*(3/4.0)
R_ANGLE_TO_WHITE_POINT = numpy.pi*(0.9653993227530295)

G_COPUNCTAL = (1.535, -.535)
G_BLIND_ANGLE = scipy.pi*(5/8.0)
G_ANGLE_TO_WHITE_POINT = numpy.pi*(0.8107602719455914)

B_COPUNCTAL = (.174, 0)
B_BLIND_ANGLE = scipy.pi*(7/8.0)
B_ANGLE_TO_WHITE_POINT = numpy.pi*(0.3734310792751017)

xyY_WHITE_POINT = (1.0/3.0, 1.0/3.0)

@click.command()
@click.option('-t', '--type', 'color_blind_type', type=click.Choice(['protanopia', 'deuteranopia', 'tritanopia']), help='Color blindness type')
@click.option('--show/--no-show', 'show_flag', is_flag=True, flag_value=True, help='Show the resulting image at the end.')
@click.option('-y', '--yes', 'yes_flag', is_flag=True, flag_value=True, help='Automatically confirm prompts.')
@click.option('--sensitivity', 'sensitivity', default=0, type=click.FLOAT, help='Color blindness sensitivity.\n0: no response\n 1: full response')

@click.option('-o', '--out', 'output_file_name', type=click.Path(exists=False, file_okay=True, dir_okay=False, resolve_path=False, writable=True), help='Set output file path. If unspecified default will be used.\n "[type]_[input_file].extension')
@click.argument('input_file_name', type=click.Path(exists=True, file_okay=True, dir_okay=False, resolve_path=False) )
def simulate(color_blind_type, sensitivity, show_flag, yes_flag, output_file_name, input_file_name):
	###
	# Check output_file_name / output format
	if (output_file_name is None):
		(head, tail) = os.path.split(input_file_name)
		(name, extension) = os.path.splitext(tail)

		output_file_name = head + '/' + str(color_blind_type) + '_' + name + '.png'
	else:
		# alwasy save as .png
		(name, extension) = os.path.splitext(output_file_name)
		output_file_name = name + '.png'


	###
	# print options/arguments
	print 'Color blind type = ' + str(color_blind_type)
	sensitivity = numpy.clip(sensitivity, 0, 1)
	print 'Sensitivity = ' + str(sensitivity)
	print 'Show resulting image? = ' + str(show_flag)
	print 'Autoconfirm prompts? = ' + str(yes_flag)
	print ''
	print 'Input image = "' + str(input_file_name) + '"'
	print 'Output image = ' + str(output_file_name) + '"'
	print ''



	###
	# read image
	image = skimage.io.imread(input_file_name, as_grey=False)

	# convert to 8-bit
	image = skimage.img_as_ubyte(image, force_copy=False)

	# remove alpha channel if present
	image = image[:,:,0:3]


	###
	# convert to CIE 1931 XYZ
	image_XYZ = skimage.color.convert_colorspace(image, fromspace='RGB', tospace='XYZ')
	image_XYZ[image_XYZ <= 0] = .001


	###
	# convert to CIE 931 xyY
	image_xy = numpy.zeros((image_XYZ.shape[0], image_XYZ.shape[1], 2), dtype=float)
	image_Y = numpy.zeros(image_XYZ.shape[0:2], dtype=float)
	image_xy[:,:,0] = image_XYZ[:,:,0] / (image_XYZ[:,:,0] + image_XYZ[:,:,1] + image_XYZ[:,:,2])
	image_xy[:,:,1] = image_XYZ[:,:,1] / (image_XYZ[:,:,0] + image_XYZ[:,:,1] + image_XYZ[:,:,2])
	image_Y = image_XYZ[:,:,1]


	###
	# processing
	mod_image_xy = numpy.zeros(image_xy.shape, dtype=float)
	mod_image_Y = numpy.zeros(image_Y.shape, dtype=float)

	if (color_blind_type == 'protanopia'): 
		copunctal = numpy.array(R_COPUNCTAL, dtype=float)
		blind_color_angle = R_BLIND_ANGLE
		angle_to_white_point = R_ANGLE_TO_WHITE_POINT

	elif (color_blind_type == 'deuteranopia'):
		copunctal = numpy.array(G_COPUNCTAL, dtype=float)
		blind_color_angle = G_BLIND_ANGLE
		angle_to_white_point = G_ANGLE_TO_WHITE_POINT

	elif (color_blind_type == 'tritanopia'):
		copunctal = numpy.array(B_COPUNCTAL, dtype=float)
		blind_color_angle = B_BLIND_ANGLE
		angle_to_white_point = B_ANGLE_TO_WHITE_POINT

	else:
		print 'Invalid color_blind_type: ' + str(color_blind_type)
		exit(1)

	# for each pixel/color value. (in xy chromatic space)
	for i,j in numpy.ndindex(image_xy.shape[0:2]):
		# find displacement/distance from copunctal poitn
		disp_from_copunctal = image_xy[i,j] - copunctal
		dist_from_copunctal = numpy.linalg.norm(disp_from_copunctal)

		# find closest point to white point along the line between the copunctal and pixelvalue.
		# found via projection
		closest_disp_from_copunctal = numpy.vdot(xyY_WHITE_POINT - copunctal, disp_from_copunctal) * (disp_from_copunctal / numpy.square(dist_from_copunctal))
		closest_dist_from_copunctal = numpy.linalg.norm(closest_disp_from_copunctal)
		closest = copunctal + closest_disp_from_copunctal

		# find displacement/distance from closest white point
		disp_from_closest = image_xy[i,j] - closest
		dist_from_closest = numpy.linalg.norm(disp_from_closest)

		# check if on blind side
		on_blind_side = onBlindSide(image_xy[i,j], color_blind_type)		


		# rescale colors based on our <senstivity> and <on_blind_side>
		# and calculate the new color values.
		simdalton_value = SimDaltonMapping(image_xy[i,j], color_blind_type)

		if (on_blind_side):
			new_dist = sensitivity * dist_from_closest
			new_color_value = closest + new_dist*(disp_from_closest/dist_from_closest)
			#new_color_value = simdalton_value*(1-sensitivity) + image_xy[i,j]*(sensitivity)
		else:
			new_color_value = simdalton_value*(1-sensitivity) + image_xy[i,j]*(sensitivity)


		# store new color value in modified xyY image.
		mod_image_xy[i,j] = new_color_value

		## DEBUG info
		# if (i== 40 and j==60) or (i==40 and j==380) or (i==40 and j==225):
		# 	print '+'
		# 	print image_xy[i,j]
		# 	print 'copunc', copunctal
		# 	print 'disp', disp_from_copunctal
		# 	print 'distance', dist_from_copunctal
		# 	print 'abs angle', abs_angle
		# 	print 'blind angle', blind_color_angle
		# 	print 'angle diff', angle_diff
		# 	print 'closest', closest
		# 	#print 'non-blind sens', non_blind_sensitivity
		# 	print 'new val', new_color_value
		# 	print 'simdalton val', simdalton_value
		# 	print '+'

	# We do not modify the luminances
	mod_image_Y = image_Y


	###
	# convert back to CIE 1931 XYZ
	mod_image_XYZ = numpy.zeros(image_XYZ.shape, dtype=float)
	mod_image_XYZ[:,:,0] = (mod_image_Y/mod_image_xy[:,:,1]) * mod_image_xy[:,:,0] # (Y/y) * x
	mod_image_XYZ[:,:,1] = mod_image_Y
	mod_image_XYZ[:,:,2] = (mod_image_Y/mod_image_xy[:,:,1]) * (1 - mod_image_xy[:,:,0] - mod_image_xy[:,:,1]) # (Y/y) * (1 -x -y)


	###
	# convert back to RGB
	mod_image = skimage.color.convert_colorspace(mod_image_XYZ, fromspace='XYZ', tospace='RGB')
	mod_image = skimage.img_as_ubyte(mod_image)


	###
	# Save to file
	save = False
	if (os.path.exists(output_file_name)):
		print 'Output file "' + str(output_file_name) + '" already exists. '

		if (yes_flag == True):
			print 'Autoconfirming overwrite.'
			save = True
		else:
			if (click.confirm('Overwrite?')):
				save = True
			else:
				print 'Aborting write to file.'
	else:
		save = True

	if (save == True):
		skimage.io.imsave(output_file_name, mod_image )
		print 'Modified image writen to = "' + str(output_file_name) + '"'

	print ''




	###
	# display results
	if (show_flag == True):
		pyplot.figure(0)
		skimage.io.imshow(image)
		pyplot.title('original: '+ str(input_file_name))
		pyplot.figure(1)
		skimage.io.imshow(mod_image)
		pyplot.title('modified: ' + str(color_blind_type) + ', sensitivity=' + str(sensitivity) + ', ' + str(output_file_name))
		pyplot.show()



# Sim Dalton method of mapping xy colors
# Proto mapping 
# start point =		(0.115807, 0.073581)
# end point = 		(0.471899, 0.527051)
#
# Deuter mapping
# start point = 	(0.102776, 0.102864)
# end point = 		(0.505845, 0.493211)
#
# Trita mapping
# start point = 	(0.045391, 0.294976)
# end point = 		(0.665764, 0.334011)
#
# Info on simdalton method:
# http://colorlab.wickline.org/colorblind/colorlab/engine.js
# *apparently the perceived color lines for the 3 dichromatic types were inferred (and linearized)
# from emprical data based on color-blind tests/test-subjects.

R_SIMDALTON_START_POINT = numpy.array([0.115807, 0.073581], dtype=float)
R_SIMDALTON_END_POINT = numpy.array([0.471899, 0.527051], dtype=float)
R_SIMDALTON_SLOPE = ( (R_SIMDALTON_END_POINT - R_SIMDALTON_START_POINT)[1] ) / ( (R_SIMDALTON_END_POINT - R_SIMDALTON_START_POINT)[0] )
R_SIMDALTON_YINT =  R_SIMDALTON_START_POINT[1] - (R_SIMDALTON_SLOPE * R_SIMDALTON_START_POINT[0])

G_SIMDALTON_START_POINT = numpy.array([0.102776, 0.102864], dtype=float)
G_SIMDALTON_END_POINT = numpy.array([0.505845, 0.493211], dtype=float)
G_SIMDALTON_SLOPE = ( (G_SIMDALTON_END_POINT - G_SIMDALTON_START_POINT)[1] ) / ( (G_SIMDALTON_END_POINT - G_SIMDALTON_START_POINT)[0] )
G_SIMDALTON_YINT = G_SIMDALTON_START_POINT[1] - (G_SIMDALTON_SLOPE * G_SIMDALTON_START_POINT[0])

B_SIMDALTON_START_POINT = numpy.array([0.045391, 0.294976], dtype=float)
B_SIMDALTON_END_POINT = numpy.array([0.665764, 0.334011], dtype=float)
B_SIMDALTON_SLOPE = ( (B_SIMDALTON_END_POINT - B_SIMDALTON_START_POINT)[1] ) / ( (B_SIMDALTON_END_POINT - B_SIMDALTON_START_POINT)[0] )
B_SIMDALTON_YINT = B_SIMDALTON_START_POINT[1] - (B_SIMDALTON_SLOPE * B_SIMDALTON_START_POINT[0])

# xy_vector := numpy array of shape (2,). Assumes dtype=float.
# color_blind_type := a string specifying color blind type
def SimDaltonMapping(xy_vector, color_blind_type):
	# check color blind type
	if (color_blind_type == 'protanopia'):
		copunctal = R_COPUNCTAL
		simdalton_slope = R_SIMDALTON_SLOPE
		simdalton_yint = R_SIMDALTON_YINT

	elif (color_blind_type == 'deuteranopia'):
		copunctal = G_COPUNCTAL
		simdalton_slope = G_SIMDALTON_SLOPE
		simdalton_yint = G_SIMDALTON_YINT

	elif (color_blind_type == 'tritanopia'):
		copunctal = B_COPUNCTAL
		simdalton_slope = B_SIMDALTON_SLOPE
		simdalton_yint = B_SIMDALTON_YINT

	else:
		print 'Invalid color_blind_type: ' + str(color_blind_type)
		exit(1)

	# compute slope and y-int of confusion line
	disp_vector = xy_vector - copunctal
	confusion_line_slope = disp_vector[1]/disp_vector[0]
	confusion_line_yint = xy_vector[1] - (confusion_line_slope * xy_vector[0])

	# compute intersecttion point
	x = (simdalton_yint - confusion_line_yint) / (confusion_line_slope - simdalton_slope)
	y = (confusion_line_slope * x) + confusion_line_yint

	return numpy.array([x,y], dtype=float)

def onBlindSide(xy_vector, color_blind_type):
	if (color_blind_type == 'protanopia'): 
		copunctal = numpy.array(R_COPUNCTAL, dtype=float)
		blind_color_angle = R_BLIND_ANGLE
		angle_to_white_point = R_ANGLE_TO_WHITE_POINT

	elif (color_blind_type == 'deuteranopia'):
		copunctal = numpy.array(G_COPUNCTAL, dtype=float)
		blind_color_angle = G_BLIND_ANGLE
		angle_to_white_point = G_ANGLE_TO_WHITE_POINT

	elif (color_blind_type == 'tritanopia'):
		copunctal = numpy.array(B_COPUNCTAL, dtype=float)
		blind_color_angle = B_BLIND_ANGLE
		angle_to_white_point = B_ANGLE_TO_WHITE_POINT

	else:
		print 'Invalid color_blind_type: ' + str(color_blind_type)
		exit(1)


	disp_from_copunctal = xy_vector - copunctal

	# find angle of the confusion line.
	# check if line is vertical
	if (disp_from_copunctal[0] == 0):
		angle = (numpy.pi/2) if (disp_from_copunctal[1]>0) else (numpy.pi*(3/2))
	else:
		angle = numpy.arctan(disp_from_copunctal[1]/disp_from_copunctal[0])
	# account for left hemisphere of circle since arctan() output is only defined from [-pi/2, pi/2]
	if(disp_from_copunctal[0] < 0):
		angle = scipy.pi + angle
	else:
		angle = angle # do nothing

	# abs_angle is restricted to be from [0, 2pi]
	abs_angle = numpy.fmod(angle + 2*numpy.pi, 2*numpy.pi)
	
	# compare confusion line angle to angle towards white point.
	angle_diff = abs_angle - angle_to_white_point

	# depending on our color blind type, and our angle relative to the white point angle
	# find out if we are on the "blind" side of the xyY colorspace
	if(color_blind_type == 'protanopia' or color_blind_type == 'deuteranopia'):
		on_blind_side = (True) if (angle_diff<=0) else (False)
	elif (color_blind_type == 'tritanopia'):
		on_blind_side = (True) if (angle_diff>=0) else (False)
	else:
		print 'Invalid color_blind_type: ' + str(color_blind_type)
		exit(1)

	return on_blind_side


# MAIN
if __name__ == '__main__':
	simulate()
