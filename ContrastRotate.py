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

##
from SimulateColorBlind import SimDaltonMapping
from SimulateColorBlind import onBlindSide


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
def contrast_rotate(color_blind_type, sensitivity, show_flag, yes_flag, output_file_name, input_file_name):
	###
	# Check output_file_name / output format
	if (output_file_name is None):
		(head, tail) = os.path.split(input_file_name)
		(name, extension) = os.path.splitext(tail)

		output_file_name = head+ '/correct_' + str(color_blind_type)[0:6] + '_' + name + '.png'
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
	# if sensitivyt is 0. we cannot perform correction
	if sensitivity == 0:
		print 'Sensitivity == 0, cannot correct color blindness'
		exit(1)


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
		simdalton_value = SimDaltonMapping(image_xy[i,j], color_blind_type)

		# calculate how much to rotate and stretch the color value
		# based on sensitivity value
		# TODO: currently uses a stretch "strength" /rotate "strength" of .3
		# Find optimal values to use depending on color blind type and/or whether the current color is on the blind side or not
		# (since on the blind side, a rotate is better than the stretch.)
		stretch = .3 * (1-sensitivity)
		rotate = .3 * (1-sensitivity)


		# compute how much of the rotation/stretch contributes to the final color value
		# based on how close to the simdalton value our original color is.
		dist_from_simdalton = numpy.linalg.norm(image_xy[i,j] - simdalton_value)
		rotate_weight = (2/numpy.pi)*numpy.arctan(2*dist_from_simdalton)
		stretch_weight = 1 - rotate_weight

		# compute new color value
		# add stretch component
		disp_from_simdalton = image_xy[i,j] - simdalton_value
		stretched_color = image_xy[i,j] + (stretch*stretch_weight)*(disp_from_simdalton/dist_from_simdalton)

		# add rotate component
		# should be pos (rotate CCW) if prota/deuter
		# should be neg (rotate CW) if tritanopia
		added_angle = rotate/(2*numpy.pi*dist_from_simdalton) 
		if (color_blind_type == 'tritanopia'):
			added_angle = -1 * added_angle
		
		# find current angle
		# check for vectical angle
		if (disp_from_simdalton[0] == 0):
			angle = (numpy.pi/2) if (disp_from_simdalton[1]>0) else (numpy.pi*(3/2))
		else:
			angle = numpy.arctan(disp_from_simdalton[1]/disp_from_simdalton[0])
		# account for left hemisphere of circle since arctan() output is only defined from [-pi/2, pi/2]
		if(disp_from_simdalton[0] < 0):
			angle = scipy.pi + angle
		else:
			angle = angle # do nothing
		# restrict angle to be from [0, 2pi]
		angle = numpy.fmod(angle + 2*numpy.pi, 2*numpy.pi)

		stretched_color_dist_from_simdalton = numpy.linalg.norm(stretched_color - simdalton_value)
		result_angle = angle + added_angle
		result_x_from_simdalton = stretched_color_dist_from_simdalton * numpy.cos(result_angle)
		result_y_from_simdalton = stretched_color_dist_from_simdalton * numpy.sin(result_angle)

		new_color_value = numpy.zeros([2], dtype=float)
		new_color_value[0] = simdalton_value[0] + result_x_from_simdalton
		new_color_value[1] = simdalton_value[0] + result_y_from_simdalton

		# store new color value in modified xyY image.
		dist_from_white = numpy.linalg.norm(image_xy[i,j] - xyY_WHITE_POINT)
		if (onBlindSide(image_xy[i,j], color_blind_type) and dist_from_white >.03):
			mod_image_xy[i,j] = new_color_value
		else:
			mod_image_xy[i,j] = image_xy[i,j]


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




# MAIN
if __name__ == '__main__':
	contrast_rotate()
