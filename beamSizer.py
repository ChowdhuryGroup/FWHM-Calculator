from heapq import nsmallest
import numpy as np
import csv

def fwhm(y_values_temp, x_values):
    y_values, temp_l, temp_r = [], [], []

    # To make 'y_values_temp', a numpy array, into a python list
    for x in range(0,len(y_values_temp)):
        y_values.append(y_values_temp[x])
    peak_height = max(y_values)
    half_peak_height = max(y_values)/2

    
    # Splitting the y_values data into before and after x_value at peak height
    y_l_temp = y_values[0:y_values.index(peak_height)]
    y_r_temp = y_values[y_values.index(peak_height):len(y_values)]

    # Finds 1st closest value to half_peak_height in y_l and y_r
    y_l = nsmallest(1, y_l_temp, key=lambda x: abs(x-half_peak_height))
    y_r = nsmallest(1, y_r_temp, key=lambda x: abs(x-half_peak_height))

    # Gets x_value pairs for y_l and y_r
    temp_l.append(x_values[y_l_temp.index(y_l[0])])

    temp_r.append(x_values[y_r_temp.index(y_r[0]) + len(y_l_temp) -1])
    fwhm_n = temp_l[0] - temp_r[0]
    return abs(fwhm_n)

xDat = []
yDat = []
with open("fwhm.tsv") as fd:
    rd = csv.reader(fd, delimiter="\t", quotechar='"')
    for row in rd:
        xDat.append(float(row[0]))
        yDat.append(float(row[1]))
print(fwhm(yDat,xDat))
