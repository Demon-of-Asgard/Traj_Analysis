#!/home/demon/anaconda3/bin/python
import os
import matplotlib.pyplot as plt
import numpy as np
def plot():
	data_folder = '../OutFiles/Traj_with_xnu/'
	for file in os.listdir(data_folder):
		#print(file)
		if(file.endswith(".datx")):
			fig = plt.figure()
			data = np.genfromtxt(data_folder+file, delimiter='\t')
			plt.plot(abs(data[:,1])/100., data[:,2], "C0-", label = "$\\nu_x$")
			plt.plot(abs(data[:,1])/100., data[:,3], "C1-", label = "$\\nu_e$")
			plt.plot(abs(data[:,1])/100., data[:,4], "C2-", label = "$\\bar\\nu_e$")
			plt.plot(abs(data[:,1])/100., data[:,5], "C1-.", label = "$bin-\\nu_e$")
			plt.plot(abs(data[:,1])/100., data[:,6], "C2-.", label = "$bin-\\bar\\nu_e$")
			plt.legend(loc=0, frameon='True')
			plt.xlabel("abs(z)[100 km]")
			plt.ylabel("$density~~[cm^{-3}]$")
			plt.title("invr2_lin-lin-"+file[:-5])
			plt.savefig(data_folder+"/plots/lin_inverser2_interpolation/"+file[:-5]+"loglin.png")


	for file in os.listdir(data_folder):
		if(file.endswith(".datx")):
			fig = plt.figure()
			plt.yscale("log")
			data = np.genfromtxt(data_folder+file, delimiter='\t')
			plt.plot(abs(data[:,1])/100., data[:,2], "C0-", label = "$\\nu_x$")
			plt.plot(abs(data[:,1])/100., data[:,3], "C1-", label = "$\\nu_e$")
			plt.plot(abs(data[:,1])/100., data[:,4], "C2-", label = "$\\bar\\nu_e$")
			plt.plot(abs(data[:,1])/100., data[:,5], "C1-.", label = "$bin-\\nu_e$")
			plt.plot(abs(data[:,1])/100., data[:,6], "C2-.",label = "$bin-\\bar\\nu_e$")
			plt.legend(loc=0, frameon='True')
			plt.xlabel("abs(z)[100 km]")
			plt.ylabel("$density~~[cm^{-3}]$")
			plt.title("invr2_lin-log-" + file[:-5])
			plt.savefig(data_folder+"/plots/lin_inverser2_interpolation/"+file[:-5]+"loglinlog.png")


def main():
	plot()


if __name__ == "__main__":
	main()
