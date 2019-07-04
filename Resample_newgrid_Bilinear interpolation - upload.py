# -*- coding: utf-8 -*-
import arcpy as ap
import numpy as np
from arcpy import env
import time
from multiprocessing.dummy import Pool as ThreadPool


env.workspace = "{workspace}"
raster_CHM = "{file}"
NetPoint = "{file}"
ra2array = ap.RasterToNumPyArray(raster_CHM)
start_time = time.time()
#==================================##==================================#
Cv = []
pool = ThreadPool(4)
#==================================##==================================#
desc = ap.Describe(raster_CHM)
height = desc.height
width = desc.width
MCH = desc.meanCellHeight
MCW = desc.meanCellWidth
#==================================##==================================#
top = ap.GetRasterProperties_management(raster_CHM, "TOP")
left = ap.GetRasterProperties_management(raster_CHM, "LEFT")
right = ap.GetRasterProperties_management(raster_CHM, "RIGHT")
bottom = ap.GetRasterProperties_management(raster_CHM, "BOTTOM")
dt = np.dtype([('extent', np.float64)])
a = np.array([top, left, right, bottom], dtype = dt)
#==================================##==================================#
cursor = ap.da.SearchCursor(NetPoint, ["SHAPE@XY"])
list_x = []
list_y = []
for row in cursor:
    x = row[0][0]
    y = row[0][1]
    list_x = np.append(list_x, x)
    list_y = np.append(list_y, y)
# ==================================##==================================#

#+++++++++++++++++++++++++++++start#==================================#
def Location_DP_model(index):
    x_p = list_x[index]
    y_p = list_y[index]
    #print(x_p, y_p)
    delta_x = x_p - a[1][0][0]
    delta_y = y_p - a[3][0][0]
    Mutiple_x = delta_x / MCW
    Mutiple_y = delta_y / MCH
    Decimal_x = abs(Mutiple_x - (int(Mutiple_x)))
    Decimal_y = abs(Mutiple_y - (int(Mutiple_y)))
    if Decimal_x <= 0.5 and Decimal_y <= 0.5:
        one = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 1)]
        two = [(int(Mutiple_x) + 0), (int(Mutiple_y) + 1)]
        three = [(int(Mutiple_x) + 0), (int(Mutiple_y) + 0)]
        four = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 0)]
    elif Decimal_x <= 0.5 and Decimal_y >= 0.5:
        one = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 1)]
        two = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 2)]
        three = [(int(Mutiple_x) + 0), (int(Mutiple_y) + 2)]
        four = [(int(Mutiple_x) + 0), (int(Mutiple_y) + 1)]
    elif Decimal_x >= 0.5 and Decimal_y >= 0.5:
        one = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 1)]
        two = [(int(Mutiple_x) + 2), (int(Mutiple_y) + 1)]
        three = [(int(Mutiple_x) + 2), (int(Mutiple_y) + 2)]
        four = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 2)]
    elif Decimal_x >= 0.5 and Decimal_y <= 0.5:
        one = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 1)]
        two = [(int(Mutiple_x) + 1), (int(Mutiple_y) + 0)]
        three = [(int(Mutiple_x) + 2), (int(Mutiple_y) + 0)]
        four = [(int(Mutiple_x) + 2), (int(Mutiple_y) + 1)]
    #===================
    return (one, two, three, four)

def GCV(FourPlocation):  #这里可能还要修改
    V_one_x = FourPlocation[0][0]
    V_one_y = (height) - FourPlocation[0][1]
    V_one = ra2array[V_one_y][V_one_x]
    V_two_x = FourPlocation[1][0]
    V_two_y = (height) - FourPlocation[1][1]
    V_two = ra2array[V_two_y][V_two_x]
    V_three_x = FourPlocation[2][0]
    V_three_y = (height) - FourPlocation[2][1]
    V_three = ra2array[V_three_y][V_three_x]
    V_four_x = FourPlocation[3][0]
    V_four_y = (height) - FourPlocation[3][1]
    V_four = ra2array[V_four_y][V_four_x]
    return (V_one, V_two, V_three, V_four)

def Bilinear_interpolation(GCV):
    x_p = list_x[0]
    y_p = list_y[0]
    # print(x_p, y_p)
    delta_x = x_p - a[1][0][0]
    delta_y = y_p - a[3][0][0]
    Mutiple_x = delta_x / MCW
    Mutiple_y = delta_y / MCH
    Decimal_x = abs(Mutiple_x - (int(Mutiple_x)))
    Decimal_y = abs(Mutiple_y - (int(Mutiple_y)))
    Decimal_x_t = Decimal_x * 0.5
    Decimal_y_t = Decimal_y * 0.5
    if Decimal_x >= 0.5 and Decimal_y >= 0.5:
        alpha = abs((Decimal_x_t - 0.25)) / 0.5
        beta = abs((Decimal_y_t - 0.25)) / 0.5
    elif Decimal_x <= 0.5 and Decimal_y <=0.5:
        alpha = abs((0.25 - Decimal_x_t)) / 0.5
        beta = abs((0.25 - Decimal_y_t)) / 0.5
    elif Decimal_x <= 0.5 and Decimal_y >=0.5:
        alpha = abs((Decimal_y_t - 0.25)) / 0.5
        beta = abs((0.25 - Decimal_x_t)) / 0.5
    elif Decimal_x >= 0.5 and Decimal_y <=0.5:
        alpha = abs((0.25 - Decimal_y_t)) / 0.5
        beta = abs((Decimal_x_t - 0.25)) / 0.5
    g_u_v = ((1 - alpha) * (1 - beta) * GCV[0]) + (alpha * (1 - beta) * GCV[1]) + \
            (beta * (1 - alpha) * GCV[2]) + (alpha * beta * GCV[3])
    print(g_u_v)
    return (g_u_v)


def workflow(index):
    FourPlocation = Location_DP_model(index)
    GCV_4 = GCV(FourPlocation)
    GUV = Bilinear_interpolation(GCV_4)
    return (GUV)
map_pro = pool.map(workflow, range(len(list_x)))

Cv = list(map_pro)
list_all = np.vstack((list_x, list_y, Cv)).transpose()
np.savetxt("Out_file", list_all, delimiter=',')