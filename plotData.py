#!/usr/local/Cellar/python3

import struct
from numpy import *
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib

experiment = "1c"

rsamp = 2000 # Hz
tstop = 90000 # ms
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

# 1c. cli and ecl vs. t for different mkcc2 values
def plot_1c():
	t = [ j/2 for j in range(0, int(tstop*2 + 1), int(rstep*2)) ]
	for i in range(0, 5):
		mkcc2_value = 0.05 + 0.1*i

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
		for l in range(8, 8 + num_samples*num_bytes, num_bytes):
			file.seek(l)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			cli_values.append(point)

		for m in range(8 + num_samples*num_bytes, 8 + 2*num_samples*num_bytes, num_bytes):
			file.seek(m)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			ecl_values.append(point)
		
		# Print the final cli and ecl values for each value of mkcc2
		print("mkcc2: %.2f" % mkcc2_value, "final cli: %.7f" % cli_values[-1], " final ecl: %.6f" % ecl_values[-1])

		plt.figure(1, figsize=(15, 10))
		plt.subplot(2, 1, 1)
		plt.plot(t[0:10000], cli_values[0:10000], '-', label="mkcc2 = %g" % mkcc2_value, lw=0.5)
		plt.subplot(2, 1, 2)
		plt.ylim([-100, -60])
		plt.plot(t[0:10000], ecl_values[0:10000], '-', label="mkcc2 = %g" % mkcc2_value, lw=0.5)

	matplotlib.rcParams.update({'font.size': 10})
	plt.subplot(2, 1, 1)
	plt.legend()
	plt.xlabel("Time (s)")
	plt.ylabel("Internal Cl- Concentration (mM)")
	plt.subplot(2, 1, 2)
	plt.legend()
	plt.xlabel("Time (s)")
	plt.ylabel("Cl- Reversal Potential (mV)")
	plt.savefig("data/1c/data.png")

# 2a. voltage vs. t for some value of egaba, for pre-post and post-pre
def plot_2a():
	t = [ j/2 for j in range(0, int(tstop*2), int(rstep*2)) ]
	cli_value = 8.00

	# Open the files for reading
	filename_pre = "data/2a/data_%.2f_pre-post.dat" % cli_value
	filename_post = "data/2a/data_%.2f_post-pre.dat" % cli_value

	def get_v_values(filename):
		file = open(filename, 'rb')

		# Get information from the header
		header_info = get_header_info(file)
		num_samples = header_info[0]
		num_bytes = header_info[1]
		val_type = header_info[2]

		v_values = []
		for l in range(8, 8 + num_samples*num_bytes, num_bytes):
			file.seek(l)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			v_values.append(point)
		return v_values

	v_values_pre = get_v_values(filename_pre)
	v_values_post = get_v_values(filename_post)

	matplotlib.rcParams.update({'font.size': 10})
	plt.figure(1, figsize=(15, 10))
	plt.subplot(2, 1, 1)
	plt.plot(t[39800:40600], v_values_pre[39800:40600], '-r', label="cli = %g" % cli_value, lw=0.5)
	plt.xlabel("Time (ms)")
	plt.ylabel("Voltage (mV)")
	plt.subplot(2, 1, 2)
	plt.plot(t[39800:40600], v_values_post[39800:40600], '-r', label="cli = %g" % cli_value, lw=0.5)
	plt.xlabel("Time (ms)")
	plt.ylabel("Voltage (mV)")
	plt.savefig("data/2a/data_v.png")

# 2b. cai vs. t for some value of egaba, for pre-post and post-pre
def plot_2b():
	t = [ j/2 for j in range(0, int(tstop*2), int(rstep*2)) ]
	cli_value = 8.00

	# Open the files for reading
	filename_pre = "data/2b/data_%.2f_pre-post.dat" % cli_value
	filename_post = "data/2b/data_%.2f_post-pre.dat" % cli_value

	def get_cai_values(filename):
		file = open(filename, 'rb')

		# Get information from the header
		header_info = get_header_info(file)
		num_samples = header_info[0]
		num_bytes = header_info[1]
		val_type = header_info[2]

		cai_values = []
		for l in range(8, 8 + num_samples*num_bytes, num_bytes):
			file.seek(l)
			data = file.read(num_bytes)
			point = struct.unpack(val_type, data)[0]
			cai_values.append(point)
		return cai_values

	cai_values_pre = get_cai_values(filename_pre)
	cai_values_post = get_cai_values(filename_post)

	matplotlib.rcParams.update({'font.size': 10})
	plt.figure(1, figsize=(15, 10))
	plt.subplot(2, 1, 1)
	plt.plot(t[20000:120000], cai_values_pre[20000:120000], '-r', label="cli = %g" % cli_value, lw=0.5)
	plt.xlabel("Time (ms)")
	plt.ylabel("Internal Ca2+ Concentration (mM)")
	plt.subplot(2, 1, 2)
	plt.plot(t[20000:120000], cai_values_post[20000:120000], '-r', label="cli = %g" % cli_value, lw=0.5)
	plt.xlabel("Time (ms)")
	plt.ylabel("Internal Ca2+ Concentration (mM)")
	plt.savefig("data/2b/data_cai.png")

if experiment == "1c":
	plot_1c()
elif experiment == "2a":
	plot_2a()
elif experiment == "2b":
	plot_2b()