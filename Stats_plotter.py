import numpy
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from matplotlib import dates
import datetime
from astropy.time import Time

plot_name = "all_eor_data_stats.pdf"
data = numpy.loadtxt("./statistics_output/All_EoR_Stats.txt")
fontsize = 30
ticksize = 25
#statistics = numpy.loadtxt("/home/ronniyjoseph/Downloads/mwa_health.csv", skiprows =1, delimiter = ',' )

#dipole_stats = numpy.loadtxt("../broken_dipole_count_2015_PhaseI.txt")[1:]

obsids = data[:, 0]
print(len(obsids))
amans_list = numpy.loadtxt("Ultimate-EOR-obsids-2014-19-full.txt")
print(len(amans_list))
gps_stamp = Time(data[:, 0], format = 'gps')
years = gps_stamp.decimalyear

print("Making plot")

figure, axes  = pyplot.subplots( 1, 1, figsize = (15,5), subplot_kw = dict(rasterized = True))
figure.subplots_adjust(bottom=0.2)

axes.scatter(years, data[:, 1], label = " Dipole", s= 10)#, "C0", label = "XX")
axes.scatter(years, data[:, 2], label = "2 Dipoles", s = 10)#, "C1", label = "YY")

axes.set_xlabel("Years", fontsize = fontsize)
axes.set_ylabel("Number of Tiles", fontsize = fontsize)
axes.legend(loc="upper right", fontsize = ticksize)
axes.tick_params(axis='both', which='major', labelsize=ticksize)

figure.savefig("./statistics_plots/" + plot_name)
figure.tight_layout()
