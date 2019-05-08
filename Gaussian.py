import csv as cs
from numpy import *
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.special import erf
import scipy


def formatnum(x, pos):
   return '$%.1f$x$10^{7}$' % (x/10000000)


def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2))) +\
           amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen2)**2)/((2*sigma2)**2)))


def _1gaussian(x, amp1,cen1,sigma1):
    return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp(-((x-cen1)**2)/((2*sigma1)**2)))


# Open csv file
File_name = 'Nap1 Bimodality.csv'
csvFile = open(File_name, "r")
reader = cs.reader(csvFile)
n = -1
k = 0
x, y = [], []
ax = []
rows = 5
cols = 3
gs = gridspec.GridSpec(rows, cols)
gs.update(hspace=0.6)
file_name = 'G3.pdf'
pp = PdfPages(file_name)
fig = plt.figure(figsize=(20,10))
formatter = FuncFormatter(formatnum)
for item in reader:
    if (str(item[0]) == '71' or str(item[0]) == '89' or str(item[0]) == '102' or str(item[0]) == '107' or str(item[0]) == '109' or str(item[0]) == '111') and k != 0:
        x1 = np.arange(x[0], x[-1], (x[-1] - x[0]) / 100)
        a1 = y[0]
        s = x[0]
        e = x[-1]
        a2 = y[-1]
        row = ((k-1) % 15 // cols)
        col = (k-1) % 15 % cols
        ax.append(fig.add_subplot(gs[row, col]))
        for p in range(1, len(y)):
            if y[p] >= a1:
                a1 = y[p]
            else:
                break
        for p in range(len(y) - 1, 0, -1):
            if y[p] >= a2:
                a2 = y[p]
            else:
                break
        c1 = x[y.index(a1)]
        c2 = x[y.index(a2)]
        u = x[2] - x[1]
        a, c = 0, 0
        if c1 != c2:
            x = np.array(x)
            y = np.array(y)
            ax[-1].bar(x, y, width=0.02)
            print(x,y)
            popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, x, y, p0=[a1*1.5, c1, u, a2*1.5, c2, u])
            perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
            if popt_2gauss[1] > popt_2gauss[4]:
                pars_1 = popt_2gauss[0:3]
                pars_2 = popt_2gauss[3:6]
            else:
                pars_2 = popt_2gauss[0:3]
                pars_1 = popt_2gauss[3:6]
            gauss_peak_1 = _1gaussian(x1, *pars_1)
            gauss_peak_2 = _1gaussian(x1, *pars_2)
            ax[-1].plot(x1, gauss_peak_1, "g")
            ax[-1].fill_between(x1, 0, gauss_peak_1, facecolor="green", alpha=0.5)
            ax[-1].plot(x1, gauss_peak_2, "y")
            ax[-1].fill_between(x1, 0, gauss_peak_2, facecolor="yellow", alpha=0.5)
            ax[-1].plot(x1, gauss_peak_2 + gauss_peak_1, "k")
            area1 = quad(_1gaussian, s, e, args=(pars_1[0], pars_1[1], pars_1[2]))
            area2 = quad(_1gaussian, s, e, args=(pars_2[0], pars_2[1], pars_2[2]))
            Columns = ['Peak1', 'Cen1', 'Peak2', 'Cen2', 'Ratio']
            Rows = ['Area']
            cell_text = [[int(area1[0]), int(pars_1[1]*100)/100.0, int(area2[0]), int(pars_2[1]*100)/100.0, int(area1[0] / (area2[0]+1)*100)/100.0]]
            the_table = ax[-1].table(cellText=cell_text, rowLabels=Rows, colLabels=Columns, loc='top')
            ax[-1].set_title(str(i[0]) + '-' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + ' ' + str(i[4]),
                             fontdict={'fontsize': 9}, pad=-9, loc='right')
            ax[-1].set_xlabel('m/Z')
            ax[-1].yaxis.set_major_formatter(formatter)
        else:
            c1 = (x[-1] - x[0]) / 4 + x[0]
            c2 = (x[-1] - x[0]) / 4 * 3 + x[0]
            x = np.array(x)
            y = np.array(y)
            print(x,y)
            ax[-1].bar(x, y, width=0.02)
            popt_2gauss, pcov_2gauss = curve_fit(_2gaussian, x, y, p0=[a1, c1, u, a2, c2, u])
            perr_2gauss = np.sqrt(np.diag(pcov_2gauss))
            if popt_2gauss[1] > popt_2gauss[4]:
                pars_1 = popt_2gauss[0:3]
                pars_2 = popt_2gauss[3:6]
            else:
                pars_2 = popt_2gauss[0:3]
                pars_1 = popt_2gauss[3:6]
            gauss_peak_1 = _1gaussian(x1, *pars_1)
            gauss_peak_2 = _1gaussian(x1, *pars_2)
            ax[-1].plot(x1, gauss_peak_1, "g")
            ax[-1].fill_between(x1, 0, gauss_peak_1, facecolor="green", alpha=0.5)
            ax[-1].plot(x1, gauss_peak_2, "y")
            ax[-1].fill_between(x1, 0, gauss_peak_2, facecolor="yellow", alpha=0.5)
            ax[-1].plot(x1, gauss_peak_2 + gauss_peak_1, "k")
            area1 = quad(_1gaussian, s, e, args=(pars_1[0], pars_1[1], pars_1[2]))
            area2 = quad(_1gaussian, s, e, args=(pars_2[0], pars_2[1], pars_2[2]))
            Columns = ['Peak1', 'Cen1', 'Peak2', 'Cen2', 'Ratio']
            Rows = ['Area']
            cell_text = [[int(area1[0]), int(pars_1[1]*100)/100.0, int(area2[0]), int(pars_2[1]*100)/100.0, int(area1[0] / area2[0]*100)/100.0]]
            the_table = ax[-1].table(cellText=cell_text, rowLabels=Rows, colLabels=Columns, loc='top')
            ax[-1].set_title(str(i[0]) + '-' + str(i[1]) + ' ' + str(i[2]) + ' ' + str(i[3]) + ' ' + str(i[4])+'eq',
                             fontdict={'fontsize': 9}, pad=-9, loc='right')
            ax[-1].set_xlabel('m/Z')
            ax[-1].yaxis.set_major_formatter(formatter)
        i = item
        x, y = [], []
        n = -1
        k += 1
        if k % 15 == 0:
            plt.savefig(pp, format='pdf')  # Save figure in pdf
            plt.close()  # Close the figure
            fig = plt.figure(figsize=(20, 10))  # Crate new figure
            ax = []
        if col == 0:
            ax[-1].set_ylabel('Intensity', {'fontsize': 9})
            ax[-1].yaxis.set_label_coords(-0.2, 0.5)
    elif k == 0:
        k += 1
        i = item
    if n > 0:
        if str(item[6]) != 'Unassigned':
            x.append(float(item[0]))
            y.append(float(item[3]))
    else:
        n = n + 1
plt.savefig(pp, format='pdf')
pp.close()
csvFile.close()


