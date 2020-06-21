#!/home/demon/anaconda3/bin/python

''' Creates the lines along different directions. '''
import numpy as np

r = np.linspace(50.0, 100.0, 50)
theta = 45.
time = 2.5
f = open("../test_trajectories/test_traj_" + str(theta) + ".dat", "w")

for rr in r:
	x = rr*np.sin(theta*np.pi/180.)
	y = rr*np.cos(theta*np.pi/180.)
	print('{0:.4f} {1:.4f} {2:.4f}'.format(rr, x, y))
	f.write('{:.4f}  {:.4f}  {:.4f} '.format(time, x, y))
	f.write("\n")
	
f.close()