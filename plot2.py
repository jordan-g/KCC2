#!/usr/local/Cellar/python3

import struct
from numpy import *
import scipy.interpolate
import matplotlib.pyplot as plt

ca_value = 0
values = []
kin_active_values = []
phos_active_values = []

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

	print(num_samples, precision)
	
	# Initialize a counter for the number of runs that we have read
	count = 0
	for i in range(1, 401):
		# Add current cai value to the values list
		ca_value = (10*i)*1e-7
		values.append(ca_value)
		
		# Set location in the file to start reading data for this run
		start = count*(16 + 2*num_samples*num_bytes) + 8
		
		kin_active_points = []
		phos_active_points = []
		
		# Get list of kin_active and phos_active values for given cai value
		for k in range(start, start + num_samples*num_bytes, num_bytes):
			file.seek(k)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			kin_active_points.append(point)

		for l in range(start + num_samples*num_bytes, start + 8 + 2*num_samples*num_bytes, num_bytes):
			file.seek(l)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			phos_active_points.append(point)
		
		# Add final values of kin_active and phos_active to their respective lists
		kin_active_values.append(kin_active_points[-1])
		phos_active_values.append(phos_active_points[-1])
		print("cai:", ca_value, "kin_active:", kin_active_points[-1], "phos_active:", phos_active_points[-1])
		# Increment the run counter
		count += 1

# Read corresponding data file, which populates the values, kin_active_values and phos_active_values lists
neuron_read("/Users/Jordan/Documents/CSB498/Jordan/data/2/data.dat")

# Plot and save graph
plt.plot(values, kin_active_values, 'r', label="Kinase")
plt.plot(values, phos_active_values, 'g', label="Phosphotase")
plt.title("Plot of steady state kin_active and phos_active vs. cai")
plt.xlabel("Ca2+ Concentration (mM)")
plt.ylabel("Percent Active (%)")
plt.legend()
plt.show()
plt.savefig("/Users/Jordan/Documents/CSB498/Jordan/data/2/data.png")