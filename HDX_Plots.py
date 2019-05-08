import csv as cs
from numpy import *
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors


def wood(df, State1, State2, Time_point):
    df[State1 + '_' + Time_point] = df[State1 + '_' + Time_point].astype('float')
    df[State2 + '_' + Time_point] = df[State2 + '_' + Time_point].astype('float')
    df['dif'] = df[State2 + '_' + Time_point] - df[State1 + '_' + Time_point]
    df[State1 + '_' + Time_point + '_SD'] = df[State1 + '_' + Time_point + '_SD'].astype('float')
    df[State2 + '_' + Time_point + '_SD'] = df[State2 + '_' + Time_point + '_SD'].astype('float')
    df['dif_err'] = df[State2 + '_' + Time_point + '_SD'] + df[State1 + '_' + Time_point + '_SD']
    df['dif_err'] = df['dif_err'].astype('float')
    x = []
    len = []
    for se in df['Sequence Number']:
        s = se.split('-')
        while '' in s:
            s.remove('')
        x.append((float(s[0]) + float(s[1])) / 2)
        len.append(float(s[1]) - float(s[0]))
    fig, ax = plt.subplots(figsize=(20, 10))
    ax.errorbar(x, df['dif'], xerr=len, marker='o', linestyle='', markersize=4, capsize=2)
    ax.grid(True)
    ax.axhline(0, color='black', lw=1)
    ax.set_xlabel('Sequence')
    ax.set_title('dif' + '_' + State1 + '_' + State2 + '_' + Time_point)
    plt.savefig('dif' + '_' + State1 + '_' + State2 + '_' + Time_point + '.eps', format='eps', dpi=1000)
    plt.show()
    return ax


def uptakeplot(df, Time_points1=[], States=[], cols=1, rows=1, file_name='Multi-page.pdf',
               color=['k', 'b', 'r', 'g', 'y']):
    # Crate grid for plot
    gs = gridspec.GridSpec(rows, cols)
    gs.update(hspace=0.5)
    x = []
    y = []
    yerr = []
    ax = []
    df.index = df['Sequence Number']
    pp = PdfPages(file_name)
    i = 0
    # Plot the uptake plot and save as pdf file
    fig = plt.figure(figsize=(7, 5))
    for Sequence_number in df['Sequence Number']:
        print(Sequence_number)
        n = 0
        row = (i // cols)
        col = i % cols
        print(row, col)
        ax.append(fig.add_subplot(gs[row, col]))  # Crate the subplot
        ax[-1].set_xscale("log", nonposx='clip')  # Set up log x
        ax[-1].set_ylim([0, float(df.loc[Sequence_number, 'MaxUptake'])])  # Set up y scale
        ax[-1].set_title(Sequence_number, fontdict={'fontsize': 6}, pad=-6, loc='right')  # Set title of plot
        ax[-1].tick_params(axis='both', labelsize=4, pad=1.2)
        if int(float(df.loc[Sequence_number, 'MaxUptake'])) // 5 == 0:
            stp = 1
        else:
            stp = int(float(df.loc[Sequence_number, 'MaxUptake'])) // 5
        ax[-1].set_yticklabels(list(range(0, int(float(df.loc[Sequence_number, 'MaxUptake'])) + stp * 2, stp)))
        print(list(range(0, int(float(df.loc[Sequence_number, 'MaxUptake'])), stp)))
        if row == rows - 1:
            ax[-1].set_xlabel('Time (s)', {'fontsize': 6})
        if col == 0:
            ax[-1].set_ylabel('Uptake (Da)', {'fontsize': 6})
            ax[-1].yaxis.set_label_coords(-0.2, 0.5)
        for State in States:
            n += 1
            for time in Time_points1:  # For 4 time points
                Line = State + '_' + time
                x.append(float(df.loc[Sequence_number, Line]))  # Get y number from df
                y.append(int(time))
                yerr.append(2 * float(df.loc[Sequence_number, Line + '_SD']))
            ax[-1].errorbar(y, x, yerr=yerr, marker='o', label=State, linewidth=0.7, markersize=0,
                            elinewidth=0.3, capsize=1, capthick=0.3, color=color[n - 1])
            # Plot one state on the subplot
            y = []
            x = []
            yerr = []
        if row == 0 and col == 0:
            ax[-1].legend(fontsize=4, loc='lower right', bbox_to_anchor=(0, 1.05))  # Set figure legend
        if i == cols * rows - 1:
            plt.savefig(pp, format='pdf')  # Save figure in pdf
            plt.close()  # Close the figure
            fig = plt.figure(figsize=(7, 5))  # Crate new figure
            ax = []
            i = -1
        i = i + 1
    if i == 0:
        plt.close()
    else:
        plt.savefig(pp, format='pdf')  # Save figure in pdf
        plt.close()  # Close the figure
    pp.close()  # Close the pdf file
    text = []
    return text


def heatmap(df, protien, State1, State2, Time_points, min=0., max=2.5, step=10, color="Blues", file_name='Heatmap.eps', step2=0):
    k = 0
    sec = list(df[protien])
    while np.core.numeric.NaN in sec:
        sec.remove(np.core.numeric.NaN)
    for time in Time_points:
        t1 = list(df[protien + '_' + State1 + '_' + time])
        t2 = list(df[protien + '_' + State2 + '_' + time])
        while np.core.numeric.NaN in t1:
            t1.remove(np.core.numeric.NaN)
        while np.core.numeric.NaN in t2:
            t2.remove(np.core.numeric.NaN)
        t1 = np.array(t1).astype(float)
        t2 = np.array(t2).astype(float)
        dif = t1 - t2
        if k == 0:
            t = copy(dif)
            k = k + 1
        else:
            t = np.vstack((t, dif))
    fig, ax = plt.subplots(figsize=(len(sec)*0.0612318+1.3243, 7))
    clmap = [(1.0, 1.0, 1.0)]
    if color == 'r':
        for c in range(step - 1):
            clmap.append((1.0 - (c + 1) * (1.0 / step) / 3, 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step)))
    elif color == 'g':
        for c in range(step - 1):
            clmap.append(
                (1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step) / 1.5, 1.0 - (c + 1) * (1.0 / step)))
    elif color == 'b':
        for c in range((step - 1)//2):
            clmap.append(
                (1.0 - ((255-100) / ((step - 1)//2) * (c+1))/255, 1.0 - ((255-149) / ((step - 1)//2) * (c+1))/255, 1.0 - ((255-235) / ((step - 1)//2) * (c+1))/255))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.append(
                ((100-100/step*(c+1))/255, (149-149/step*(c+1))/255,
                 1.0 - ((255 - 150) / ((step - 1)//2) * (c + 1)) / 255))
    elif color == 'rb':
        clmap = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]
        for c in range(step - 1):
            clmap.append((1.0 - (c + 1) * (1.0 / step) / 3, 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step)))
        for c in range(step2 - 1):
            clmap.insert(0,
                    (1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step) / 1.5))
    elif color == 'br':
        clmap = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]
        for c in range((step - 1)//2):
            clmap.append(
                (1.0 - ((255-100) / ((step - 1)//2) * (c+1))/255, 1.0 - ((255-149) / ((step - 1)//2) * (c+1))/255, 1.0 - ((255-235) / ((step - 1)//2) * (c+1))/255))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.append(
                ((100-100/step*(c+1))/255, (149-149/step*(c+1))/255,
                 1.0 - ((255 - 150) / ((step - 1)//2) * (c + 1)) / 255))
        for c in range(step2 - 1):
            clmap.insert(0, (1.0 - (c + 1) * (1.0 / step) / 3, 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step)))
    elif color == 'o':
        for c in range((step - 1)//2):
            clmap.append(
                (1, 1.0 - ((255-102) / ((step - 1)//2) * (c+1))/255, 1.0 - (255 / ((step - 1)//2) * (c+1))/255))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.append(((255-(230-51)/step*(c+1))/255, (102-(102-20)/step*(c+1))/255, 0))
    elif color == 'ob':
        clmap = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]
        for c in range((step - 1)//2):
            clmap.append(
                (1, 1.0 - ((255-102) / ((step - 1)//2) * (c+1))/255, 1.0 - (255 / ((step - 1)//2) * (c+1))/255))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.append(((255-(230-128)/step*(c+1))/255, (102-(102-51)/step*(c+1))/255, 0))
        for c in range(step2 - 1):
            clmap.insert(0,
                    (1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step) / 1.5))
    elif color == 'y':
        for c in range((step - 1)//2):
            clmap.append(
                (1, 1, 1.0 - (255 / ((step - 1)//2) * (c+1))/255))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.append(((255-(255-77)/step*(c+1))/255, (255-(255-77)/step*(c+1))/255, 0))
    elif color == 'gr':
        for c in range(step - 1):
            clmap.append((1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step)))
    elif color == 'bp':
        clmap = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]
        for c in range((step - 1)//2):
            clmap.append(
                (1.0 - ((255-100) / ((step - 1)//2) * (c+1))/255, 1.0 - ((255-149) / ((step - 1)//2) * (c+1))/255, 1.0 - ((255-235) / ((step - 1)//2) * (c+1))/255))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.append(
                ((100-100/step*(c+1))/255, (149-149/step*(c+1))/255,
                 1.0 - ((255 - 150) / ((step - 1)//2) * (c + 1)) / 255))
        for c in range((step - 1)//2):
            clmap.insert(0,
                (1, 1.0 - ((255-0) / ((step - 1)//2) * (c+1))/255, 1))
        for c in range(step - 1 - (step - 1) // 2):
            clmap.insert(0,
                         ((255-(225-23)/(step - 1)//2*(c+1))/255, 0, (255-(225-23)/(step - 1)//2*(c+1))/255))
    else:
        clmap = [(1.0, 1.0, 1.0), (1.0, 1.0, 1.0)]
        for c in range(step - 1):
            clmap.append((1.0 - (c + 1) * (1.0 / step) / 3, 1.0 - (c + 1) * (1.0 / step), 1.0 - (c + 1) * (1.0 / step)))
        for c in range(step2-1):
            clmap.insert(0, (75/255, 140/255, 97/255))
    cmap = mpl.colors.ListedColormap(clmap)
    im = ax.imshow(t, aspect=3, cmap=cmap, vmin=min, vmax=max)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Dif', rotation=-90, va="bottom", fontsize=5)
    cbar.set_ticks(np.linspace(min, max, step + step2 + 1))
    ax.set_xticks(np.arange(len(sec)))
    ax.set_yticks(np.arange(len(Time_points)))
    ax.set_xticklabels(sec)
    ax.set_yticklabels(Time_points)
    ax.set_facecolor('white')
    ax.tick_params(axis='x', labelsize=3.5, pad=0.9, length=3.2)
    ax.tick_params(axis='y', labelsize=10)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right", va='center', rotation_mode="anchor")
    plt.title(protien + '_' + State1 + '-' + State2, {'fontsize':4})
    fig.tight_layout()
    plt.savefig(file_name, format='eps', dpi=1000)
    plt.show()
    return k


# Open csv file
File_name = 'TRAMP HDX AT.csv'
csvFile = open(File_name, "r")
reader = cs.reader(csvFile)
# Set up time points and states
Time_points1 = ['10', '100', '1000', '10000', '100000']
States1 = ['AT comple', 'AT complex+RNA']
Proteins = ['AIR2_YEAST', 'PAP2_YEAST']
Columns = []
# Crate column to store data
for Protein in Proteins:
    Columns.append(Protein)
    for State in States1:
        Columns.append(State)
        for Time in Time_points1:
            Columns.append(Protein + '_' + State + '_' + Time)
            Columns.append(Protein + '_' + State + '_' + Time + '_SD')
Data1 = pd.DataFrame(columns=Columns)
n, i = 0, 0
Sequence = ''
Sequence_number = ''
# Read data form reader to Data
for item in reader:
    if n != 0:
        if item[0] in Proteins:
            if Sequence_number != item[1] + '-' + item[2]:
                Sequence_number = item[1] + '-' + item[2]
                Data1.loc[i, item[0]] = Sequence_number
                i += 1
            Time = str(int(float(item[9]) * 60 + 0.5))
            State = item[8]
            Protein = item[0]
            if Time != '0':
                Data1.loc[i, Protein + '_' + State + '_' + Time] = item[12]
                Data1.loc[i, Protein + '_' + State + '_' + Time + '_SD'] = item[13]
    else:
        n = n + 1
csvFile.close()
# Save Data as csv file
Data1.to_csv("For plot.csv", index=False, sep=',')
# T = uptakeplot(Data1, Time_points1, States1, 5, 4, file_name='Mtr4.pdf',
#                color=[(192 / 255, 0, 0), 'k', (192 / 255, 0, 0), (22 / 255, 54 / 255, 92 / 255),
#                       'sienna'])
# for time in Time_points1:
#     for state in States1:
#         Data1[state + '_' + time] = Data1[state + '_' + time].astype('float')
#     Data1['Sub1' + '_' + time] = Data1['Mtr4' + '_' + time] - Data1['Mtr4+RNA' + '_' + time]
#     Data1['Sub3' + '_' + time] = Data1['TRAMP Complex' + '_' + time] - Data1['TRAMP Complex+RNA' + '_' + time]
for protein in Proteins:
    K = heatmap(Data1, protein, 'AT comple', 'AT complex+RNA', Time_points1, max=5, step=10, color='r', min=0, step2=0,
            file_name='{}_Mtr4_RNA_TRAMP.eps'.format(protein))
# # for time in Time_points1:
#     H = wood(df, 'Apo', 'ADP', time)

