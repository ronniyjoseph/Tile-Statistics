import numpy
import matplotlib
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from matplotlib import dates
import datetime
from astropy.time import Time

statistics = numpy.loadtxt("../broken_tile_count_2017B_PhaseII_Compact.txt")
statistics = numpy.loadtxt("/home/ronniyjoseph/Downloads/mwa_health.csv", skiprows =1, delimiter = ',' )

dipole_stats = numpy.loadtxt("../broken_dipole_count_2015_PhaseI.txt")[1:]
fontsize = 15
#reshape = dipole_stats.reshape(4,4)
#x = numpy.arange(-2,2,1)
#xx, yy = numpy.meshgrid(x,x)
#figure = pyplot.figure(figsize = (5,5))
#ax = figure.add_subplot(111, projection = '3d')
#ax.bar3d(xx.flatten(),yy.flatten(),numpy.zeros(len(dipole_stats)),1, 1, dipole_stats, cmap = cm.get_cmap('jet'))
#ax.imshow(reshape)
#pyplot.show()

d = statistics[:, 0 ]
s = d/1000
ms = d-1000*s
dts = datetime.datetime.fromtimestamp(s[0])
#fds = dates.date2num(dts) # converted

hfmt = dates.DateFormatter('%m/%d %H:%M')
t = Time(d, format='gps')

year = t.decimalyear
print(year)
print("Making plot")
figure = pyplot.figure(figsize = (25,5))
figure.subplots_adjust(bottom=0.2)

figure.tight_layout()
axes = figure.add_subplot(111)

axes.scatter(year, statistics[:, 1], label = "XX")#, "C0", label = "XX")
axes.scatter(year, statistics[:, 3], label = "YY")#, "C1", label = "YY")

axes.set_xlabel("Observation ID number", fontsize = fontsize)
axes.set_ylabel("Number of Tiles", fontsize = fontsize)
#axes.set_title("2017", fontsize = fontsize)

#axes.tick_params(labelrotation = 20)
#axes.ticklabel_format(axis = 'x', style = 'plain')
#axes.xaxis.set_major_locator(MaxNLocator(integer=True))
#axes.set_xticks(ticks = d[::25] )
#axes.set_xticklabels(labels = d[::25] )

axes.legend(loc="upper right")

figure.savefig("simple_mwa_stats_2017.pdf")
pyplot.show()