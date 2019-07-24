import numpy
import matplotlib
from matplotlib import pyplot
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from matplotlib import dates
import datetime
from astropy.time import Time
import os

stats_folder = "/data/rjoseph/Hybrid_Calibration/Dipole_Statistics"

fontsize = 15

stats_file_list = os.listdir(stats_folder)

file_counter = 0
for file in stats_file_list:
    new_data = numpy.loadtxt(stats_folder + "/" + file)
    if file_counter == 0:
        data = new_data
    else:
        data = numpy.append(data, new_data, axis = 0)
    file_counter += 1

print(f"We have data for {len(data[:, 0])} obsids")

obsid = Time(data[:, 0], format = 'gps')
year = obsid.decimalyear
#pyplot.scatter(data[:, 0], data[:, 1])
#pyplot.show()
figure, axes = pyplot.subplots(1,1, figsize = (25,5))
axes.plot(year, data[:, 1], label = "1 Dipole", s = 10, alpha = 0.3)
axes.plot(year, data[:, 2], label = "2 Dipoles",s=10,   alpha = 0.3)

axes.set_xlabel("Year", fontsize = fontsize)
axes.set_ylabel("Number of Tiles", fontsize = fontsize)
axes.legend(loc="upper right")
#figure.savefig(stats_folder + "all_years.pdf")
pyplot.show()