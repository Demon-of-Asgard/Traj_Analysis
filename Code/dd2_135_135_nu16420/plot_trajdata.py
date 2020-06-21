#!/home/demon/anaconda3/bin/python3

#========================================================================================
import sys
import os
import matplotlib.pyplot as plt
import numpy as np

#========================================================================================

def plot(file_name, xcol, zcol, xlabel = None, ylabel = None):

    file = os.getcwd()+"/" + file_name

    trajdata = np.genfromtxt(file, skip_header =1)

    plt.plot(trajdata[:,xcol],trajdata[:,zcol],"C0-",label = "$ray-n_{\\nu_e}$")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(file)
    plt.show()

#========================================================================================

def mplot():
    try:
        file_name = sys.argv[1]
        try:
            xcol = int(sys.argv[2])
        except:
            print("x-column set to default = 0 \n")
            xcol = 0;

        try:
            ycol = int(sys.argv[3])
        except:
            print("y-column set to default = 1 \n")
            ycol = 1

    except:
        print("First is the filename. It is mandatory. Try ploting again with a valid filename. \n")
        exit()


    plot(file_name, xcol, ycol)

if __name__ == "__main__":
    mplot()


#========================================================================================
