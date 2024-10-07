import argparse
import os
from os import path
import subprocess
from time import time
# import matplotlib.pyplot as plt

SAMP_SIZES = [ 1000000, 1500000, 2500000, 278197, 3000000, 5000000, 10000000, 30000000, 50000000, 100000000 ]
NUM_SAMPLES = [ 10, 5, 5, 5, 3, 3, 3, 3, 3, 3 ]

parser = argparse.ArgumentParser(prog="", description="", epilog="")
args = None

def valid_args(args):
	result = True

	result &= path.isfile(args.old)
	result &= path.isfile(args.new)
	result &= path.isfile(args.gen)
	result &= path.isdir(args.tmp)

	if result == False:
		raise ValueError("Invalid arguments")

def parse_args():
	global args

	parser.add_argument("-e", "--exec", help="path of old fastx quality statistics", type=str, required=True)
	parser.add_argument("-g", "--gen", help="path of sample generator", type=str, required=True)
	parser.add_argument("-t", "--tmp", help="temporary sample directory", type=str, required=True)

	args = parser.parse_args()
	valid_args(args)

def gen_samp(size):
	command = [ args.gen, " -s fastq", f" --nr {size}" f" -o {path.join(args.tmp, size)}" " --mns 50", " --mxs 200" ]
	command = ''.join(command)

	subprocess.run(command, capture_output=False, text=False, check=True)

def del_samp(size):
	os.remove(path.join(args.tmp, size))

def execute_qual_stats(program_path, size):
	command = [ program_path, f" -i {path.join(args.tmp, size)}" ]
	command = ''.join(command)

	start_time = time()
	exec_res = subprocess.run(command, capture_output=True, text=True, check=True)
	end_time = time()

	return ((end_time - start_time) * 1000, exec_res.stdout)

def compare_results(old, new, size, session):
	old_output = old[1].replace('\n', '').replace('\r', '').strip() # Remove line feed & carriage return
	new_output = new[1].replace('\n', '').replace('\r', '').strip()

	print(f"[{size}:{session}] {'PASSED' if old_output == new_output else 'MISMATCH'}, Old: {old[0]}ms, New: {new[0]}ms")

def execute_comparison():
	for i in range(len(SAMP_SIZES)):
		samp_size = str(SAMP_SIZES[i])
		
		for j in range(NUM_SAMPLES[i]):
			gen_samp(samp_size)

			result = execute_qual_stats(args.old, samp_size)
			print(f"[{session}] Time: {result[0]}ms")

			del_samp(samp_size)

try:
	parse_args()
	execute_comparison()

except Exception as e:
	print(e)