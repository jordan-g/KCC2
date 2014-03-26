#!/usr/local/Cellar/python3

import struct
from numpy import *
import scipy.interpolate
import matplotlib.pyplot as plt

experiment = "1c"

rsamp = 2000 # Hz
tstop = 45000 # ms
rstep = 1/(rsamp/1000) # ms

def get_header_info(file):
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

	return [num_samples, num_bytes, val_type]

def plot_1c():
	# 1c. cli and ecl vs. t for different mkcc2 values
	t = [ j/2 for j in range(0, int(tstop*2 + 1), int(rstep*2)) ]
	print(len(t))
	for i in range(1, 11):
		mkcc2_value = 0.05*i

		# Open the file for reading
		filename = "data/1c/data_%.2f.dat" % mkcc2_value
		file = open(filename, 'rb')

		# Get information from the header
		header_info = get_header_info(file)
		num_samples = header_info[0]
		num_bytes = header_info[1]
		val_type = header_info[2]

		cli_values = []
		ecl_values = []
		
		for l in range(8, 8 + num_samples*num_bytes + 1, num_bytes):
			file.seek(l)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			cli_values.append(point)

		for m in range(8 + num_samples*num_bytes, 8 + 2*num_samples*num_bytes + 1, num_bytes):
			file.seek(m)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			ecl_values.append(point)

		plt.figure(1)
		plt.plot(t, cli_values, '-', label="mkcc2 = %g" % mkcc2_value)
		plt.figure(2)
		plt.plot(t, ecl_values, '-', label="mkcc2 = %g" % mkcc2_value)

	plt.figure(1)
	plt.savefig("data/1c/data_cli.png")
	plt.figure(2)
	plt.savefig("data/1c/data_ecl.png")

if experiment == "1c":
	plot_1c()