import shapefile
import numpy as np

file = shapefile.Reader("path.shp")
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

def read_species():
    NAME = raw_data[:, 2]
    unique_name = np.unique(NAME)
    return (unique_name)
#print(read_species())
species_name = read_species()

def rename():
    for index_rename in range(len(raw_data)):
        raw_data[index_rename][2] = np.argwhere(species_name == raw_data[index_rename][2])[0][0]
    return (raw_data)
rename_results = rename()
print(rename_results)

if __name__ == "__main__":

    np.savetxt("outpgth.csv",
               rename_results, delimiter=',', fmt = '%s')