#!/home/demon/anaconda3/bin/python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

def getinfo(file):
    stat = False
    data = np.genfromtxt(file, skip_header=1)
    thetas = data[:,10]
    for theta in thetas:
        if(theta<=45.0):
            stat = True
    return(stat)


def get_trajectories(root_PATH):
    stat = False
    files = []
    for file in os.listdir(root_PATH):
        if(file.startswith('trajectory')):
            stat = getinfo(root_PATH+file)
            if(stat==True):
                files.append(file)
    return(files)


def echo(files):
    for id,file in enumerate(files):
        print('{}-->{}'.format(id,file))


def register(files):
    np.savetxt("../OutFiles/trajectory_names.txt", np.array(files).T, fmt = ["%s"])


def plot(root_PATH, files):
    fig = plt.figure()
    for file in files:
        x = []
        z = []
        data = np.genfromtxt(root_PATH+file, skip_header=1)
        rs = data[:,9]
        thetas = data[:,10]
        for id,_ in enumerate(thetas):
            x.append(rs[id]*np.sin(thetas[id]*np.pi/180.0))
            z.append(rs[id]*np.cos(thetas[id]*np.pi/180.0))
        plt.xlim(0.0,100.0)
        plt.ylim(0.0, 100.0)
        plt.plot(x,z, '-')
    plt.xlabel('$x[km]$')
    plt.ylabel('$y[km]$')
    plt.title("$Trajectories$")
    plt.show()


def main():
    root_PATH = '../dd2_135_135_nu16420/'
    files=get_trajectories(root_PATH)
    register(files)
    #echo(files)
    plot(root_PATH, files)


if( __name__ == "__main__"):
    main()