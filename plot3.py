#!/usr/local/Cellar/python3

import struct
from numpy import *
import scipy.interpolate
import matplotlib.pyplot as plt

ca_value = 0
values = []
mkcc2_values = []

def neuron_read(filename):
	data = []
	file = open(filename, 'rb')
	
	# Get the number of data points for each run
	data = file.read(4)
	num_samples = struct.unpack('i', data)[0]

	# Get the precision argument
	file.seek(4)
	data = file.read(4)
	precision = struct.unpack('i', data)[0]
	
	# Set the right number of bytes per data point
	if precision == 1:
		num_bytes = 1
		val_type = 'c'
	elif precision == 2:
		num_bytes = 2
		val_type = 'h'
	elif precision == 3:
		num_bytes = 4
		val_type = 'f'
	elif precision == 4:
		num_bytes = 8
		val_type = 'd'
	
	# Initialize a counter for the number of runs that we have read
	count = 0
	print("---------- cai:", ca_value, "mM --------------------")
	for j in range(1, 11):
		for k in range(1, 11):
			# Add current cai, R_M and R_MP values to the values list
			values.append((ca_value, j, k))
			
			# Set location in the file to start reading data for this run
			start = count*(8 + num_samples*num_bytes) + 8
			
			mkcc2_points = []
			
			# Get list of mkcc2i values for given cai, R_M, and R_MP values
			for l in range(start, start + num_samples*num_bytes, num_bytes):
				file.seek(l)
				data = file.read(num_bytes)
				point = struct.unpack(val_type, data)[0]
				mkcc2_points.append(point)
			
			# Add final value of mkcc2i to the mkcc2_points list
			mkcc2_values.append(mkcc2_points[-1])
			print("R_M:", j, "R_MP:", k, "mkcc2:", mkcc2_points[-1])
			# Increment the run counter
			count += 1

# Loop over cai values from 1 to 10
for i in range(1, 11):
	ca_value = (50e-6)*(1 + i)
	values = []
	mkcc2_values = []
	
	# Read corresponding data file, which populates the values and mkcc2_values lists
	neuron_read("/Users/Jordan/Documents/CSB498/Jordan/data/3/data_%d.dat" % i)
	
	# set x to R_M values, y to R_MP values and z to final mkcc2 values
	x = []
	y = []
	for (l, j, k) in values:
		x.append(j)
		y.append(k)
	z = mkcc2_values
	
	# Interpolate and plot
	xi, yi = linspace(min(x), max(x), 100), linspace(min(y), max(y), 100)
	xi, yi = meshgrid(xi, yi)
	rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
	zi = rbf(xi, yi)
	plt.imshow(zi, vmin=min(z), vmax=max(z), origin='lower',
		   extent=[min(x), max(x), min(y), max(y)])
	plt.scatter(x, y, c=z)
	plt.title("Contour plot of steady state mkcc2 value vs. R_M and R_MP values for cai = %.4f" % ca_value, fontsize=10)
	plt.xlabel("R_M value")
	plt.ylabel("R_MP value")
	if i == 1:
		plt.colorbar()
	plt.savefig("/Users/Jordan/Documents/CSB498/Jordan/data/data_%d.png" % i)