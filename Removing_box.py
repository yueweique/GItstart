import shapefile
import numpy as np
from collections import Counter

file = shapefile.Reader("path.shp")
inFile_feild = np.loadtxt("Pretreated.csv",
                          dtype = np.float64, delimiter = ',', )
shapes = file.shapes()
np.set_printoptions(suppress=True)

def read_recors():
    shp_records = file.records()
    X = []
    Y = []
    NAME = []
    for index_read in range(len(shp_records)):
        X = np.append(X, shp_records[index_read][3])
        Y = np.append(Y, shp_records[index_read][4])
        NAME = np.append(NAME, shp_records[index_read][0])
    data_Arrangement = np.vstack((X, Y, NAME)).transpose()
    return (data_Arrangement)
raw_data = read_recors()
#print(raw_data)

#def creat_remoing_box():

boundary_X_MIN = file.bbox[0]
boundary_X_MAX = file.bbox[2]
boundary_Y_MIN = file.bbox[1]
boundary_Y_MAX = file.bbox[3]
#===========================================================#
Value_of_boxlength = 10
#===========================================================#
def first_P(cell_size):
    first_P_X = boundary_X_MIN + (cell_size / 2)
    first_P_Y = boundary_Y_MIN + (cell_size / 2)
    return ([first_P_X, first_P_Y])

def remove_times(cell_size):
    d = cell_size
    X_CC_min_scale = boundary_X_MIN + (d / 2)  # 定义随机立方体X主分量上的最小取值
    Y_CC_min_scale = boundary_Y_MIN + (d / 2)  # 定义随机立方体Y主分量上的最小取值
    # ===============>
    X_CC_max_scale = boundary_X_MAX - (d / 2)  # 定义随机立方体X主分量上的最大取值
    Y_CC_max_scale = boundary_Y_MAX - (d / 2)  # 定义随机立方体Y主分量上的最大取值
    # ======================================================>
    X_move_range = X_CC_max_scale - X_CC_min_scale
    Y_move_range = Y_CC_max_scale - Y_CC_min_scale
    # ======================================================>
    X_times = int(X_move_range / d) + 1
    Y_times = int(Y_move_range / d) + 1
    return (X_times, Y_times)

def remove_point(length):
    list_X_boxceter = []
    list_Y_boxceter = []

    for i_Y in range(remove_times(Value_of_boxlength)[1]):
        for i_X in range(remove_times(Value_of_boxlength)[0]):
            X_in_X = first_P(Value_of_boxlength)[0] + (length * i_X)
            Y_in_X = first_P(Value_of_boxlength)[1] + (length * i_Y)
            # ===============>
            list_X_boxceter = np.append(list_X_boxceter, X_in_X)
            list_Y_boxceter = np.append(list_Y_boxceter, Y_in_X)
            # ===============>
            list_XYZ = np.vstack((list_X_boxceter, list_Y_boxceter)).transpose()
        #print("processing: ", (((i_Y) / (remove_times(Value_of_boxlength)[1]))) * 100, '%')
    return (list_XYZ)
list_center = remove_point(Value_of_boxlength)

#print(list_center)
def workflow():
    list_results_Caculate = []
    list_results_Caculate_X = []
    list_results_Caculate_Y = []
    list_ALL_results_Caculate = []
    for index_workflow in range(len(list_center)):
        # =============clip from any cell(CAC)=========================#
        def CAC(X,Y):
            X_cube_min = X - (Value_of_boxlength / 2)
            Y_cube_min = Y - (Value_of_boxlength / 2)
            # =========================>
            X_cube_max = X + (Value_of_boxlength / 2)
            Y_cube_max = Y + (Value_of_boxlength / 2)
            #==========================>
            X_J = np.logical_and(inFile_feild[:, 0] >= X_cube_min, inFile_feild[:, 0] <= X_cube_max)
            Y_J = np.logical_and(inFile_feild[:, 1] >= Y_cube_min, inFile_feild[:, 1] <= Y_cube_max)
            XY_J = np.logical_and(X_J, Y_J)
            keep_points = inFile_feild[XY_J]
            return (keep_points)

        #test = CAC(list_center[index_workflow][0], list_center[index_workflow][1])
        #print(test)
        #print(test.dtype, test)

        def N(Total):
            N = len(Total)
            return (N)
        N_number = N(CAC(list_center[index_workflow][0], list_center[index_workflow][1]))

        def Ni(Total):
            All_SP_in_cell = Total[:, 2]
            SP_type = np.unique(All_SP_in_cell)                       #计数共有几种
            list_Ni_counts = []
            SP_Name = []
            Ni_in_cell = []
            for index_Count in range(len(SP_type)):
                mask = (All_SP_in_cell == SP_type[index_Count])
                Count_in_every_Type = All_SP_in_cell[mask]            #计数种i的个体数
                list_Ni_counts = np.append(list_Ni_counts, len(Count_in_every_Type))
                SP_Name = np.append(SP_Name, SP_type[index_Count])
                Ni_in_cell = np.vstack((SP_Name, list_Ni_counts)).transpose()
            return (Ni_in_cell)
        #print(Ni(test))
        #Ni_list = Ni(CAC(list_center[index_workflow][0], list_center[index_workflow][1]))

        # =============caculate_Biodiversity(CB)=========================#
        def CB():
            Pi_list = []
            for index_Caculate in range(len(Ni(CAC(list_center[index_workflow][0], list_center[index_workflow][1])))):
                Pi = Ni(CAC(list_center[index_workflow][0], list_center[index_workflow][1]))[index_Caculate][1] /\
                     N_number
                #print(Pi)
                Pi_1 = ((np.log10(Pi)) * Pi)
                Pi_list = np.append(Pi_list, Pi_1)
            #print(Pi_list)
            pi_sum = np.sum(Pi_list)
            H_d = - pi_sum
            return (H_d)
        #print(CB())
        every_result_in_cell = CB()
        list_results_Caculate_X = np.append(list_results_Caculate_X, list_center[index_workflow][0])
        list_results_Caculate_Y = np.append(list_results_Caculate_Y, list_center[index_workflow][1])
        list_results_Caculate = np.append(list_results_Caculate, every_result_in_cell)
        list_ALL_results_Caculate = np.vstack((list_results_Caculate_X, list_results_Caculate_Y,
                                               list_results_Caculate)).transpose()
        print("Processing: ", (((index_workflow) / (len(list_center))) * 100),'%')
    return (list_ALL_results_Caculate)
#print(workflow())
list_ALL_FIN = workflow()
np.savetxt("outpath.csv", list_ALL_FIN, delimiter=',')
