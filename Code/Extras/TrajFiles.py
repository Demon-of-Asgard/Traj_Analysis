#!/usr/bin/python3
import os
import numpy as np


def get_filenames():
	files = []
	rootpath = '../'
	for folder in os.listdir(rootpath):
		if (folder.endswith("dd2_135_135_nu16420")):
			for file in os.listdir(rootpath+folder):
				if(file.startswith("trajectory")):
					file_PATH = file 
					files.append(file_PATH)
	files = np.array(files).T
	np.savetxt("../OutFiles/trajectory_names_full.txt", files, fmt = ["%s"])
	return files

def print_filenames(filenames):
	count = 0
	for filename in filenames:
		count += 1
		print('{}---> {}'.format(count,filename))

def main():
	filenames = get_filenames()
	print_filenames(filenames)


if (__name__ == "__main__"):
	main()
