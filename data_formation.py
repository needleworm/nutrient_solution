import numpy as np
import os
import system_simulator as ss
import random
import sys

resultdir = sys.argv[1]
networkdir = sys.argv[2]
normalize = sys.argv[3]
ISE_observable = sys.argv[4]

if not resultdir.endswith("/"):
    resultdir += "/"

if not networkdir.endswith("/"):
    networkdir += "/"

reslst = os.listdir(resultdir)
modlst = os.listdir(networkdir)

# Network_files handling
networkQuery = []
for el in modlst:
    if "nutrient_plant" in el and "txt" in el and "result" not in el:
        networkQuery.append(networkdir + el)
networkQuery.sort()

print(str(len(networkQuery)) + " network query files found.")

byunghyun_coeficients = []
for el in networkQuery:
    temp_network = ss.Network(el)
    byunghyun_coeficients.append(temp_network.byunghyun_coefficients()[2:])
    # Water dissociation is not plant-ion interaction.

if ISE_observable not in "TruetrueTRUE":
    INITIAL_STATE = temp_network.Xs
    MOLECULAR_WEIGHT = temp_network.molecular_weight
else:
    i_no3 = temp_network.nameidx["[NO3-]"]
    i_ca = temp_network.nameidx["[Ca++]"]
    i_k = temp_network.nameidx["[K+]"]
    i_nh4 = temp_network.nameidx["[NH4+]"]
    i_h = temp_network.nameidx["[H+]"]
    OBSERVABLE_IDXs = [i_no3, i_ca, i_k, i_nh4, i_h]
    OBSERVABLE_IDXs.sort()
    MOLECULAR_WEIGHT = temp_network.molecular_weight[OBSERVABLE_IDXs]
    INITIAL_STATE = temp_network.Xs[OBSERVABLE_IDXs]
del(temp_network)

# Result_files handling
resultFiles = []
for el in reslst:
    if "result" in el and "csv" in el:
        resultFiles.append(resultdir + el)
resultFiles.sort()

print(str(len(resultFiles)) + " simulation result files found.")

numElements = 0
print("start data processing...\n")
DATA_tuple = []

for i, el in enumerate(resultFiles):
    txt = open(el)
    count = 0
    for line in txt:
        splt = line.strip().split(",")
        if numElements == 0:
            numElements = len(splt)
        if "[" in line:
            continue
        if len(splt) != numElements:
            print(el + " has wrong line")
            print(line)
            continue
        if ISE_observable not in "TruetrueTRUE":
            X = INITIAL_STATE - np.asarray(splt, dtype=np.float64)
        else:
            X = INITIAL_STATE - np.asarray(splt, dtype=np.float64)[OBSERVABLE_IDXs]
        TDS = X*MOLECULAR_WEIGHT
        X = np.append(X, TDS)
        Y = byunghyun_coeficients[i]
        DATA_tuple.append((X, Y))
        count += 1
    txt.close()
    print(str(count) + " data from " + el)


print("\nTotal " + str(len(DATA_tuple)) + " date prepared.\n")


random.shuffle(DATA_tuple)

num_test_data = int(len(DATA_tuple)/5)

TEST_X = []
TEST_Y = []

TRAINING_X = []
TRAINING_Y = []

for i, el in enumerate(DATA_tuple):
    X, Y = el
    if i < num_test_data:
        TEST_X.append(X)
        TEST_Y.append(Y)
    else:
        TRAINING_X.append(X)
        TRAINING_Y.append(Y)

TEST_X = np.asarray(TEST_X)
TEST_Y = np.asarray(TEST_Y)

TRAINING_X = np.asarray(TRAINING_X)
TRAINING_Y = np.asarray(TRAINING_Y)

print("The size of Training Data is : " + str(len(TRAINING_X)))
print("The size of Test Data is : " + str(len(TEST_X)))

print("DATA SAVED")

if normalize in "TRUEtrueTrue":
    print("NORMALIZATION...")
    TEST_X -= np.min(TEST_X)
    TEST_X /= np.max(TEST_X)

    TRAINING_X -= np.min(TRAINING_X)
    TRAINING_X /= np.max(TRAINING_X)

    TEST_Y = np.log10(TEST_Y)
    TEST_Y -= np.min(TEST_Y)
    TEST_Y /= np.max(TEST_Y)

    TRAINING_Y = np.log10(TRAINING_Y)
    TRAINING_Y -= np.min(TRAINING_Y)
    TRAINING_Y /= np.max(TRAINING_Y)

filename_test_x = "test_X.npy"
filename_test_y = "test_y.npy"
filename_training_x = "training_X.npy"
filename_training_y = "training_Y.npy"

if ISE_observable in "TruetrueTRUE":
    filename_test_x = "ISE_obs_" + filename_test_X
    filename_test_y = "ISE_obs_" + filename_test_y
    filename_training_x = "ISE_obs_" + filename_training_x
    filename_training_y = "ISE_obs_" + filename_training_y

if normalize in "TruetrueTRUE":
    filename_test_x = "Normalized_" + filename_test_X
    filename_test_y = "Normalized_" + filename_test_y
    filename_training_x = "Normalized_" + filename_training_x
    filename_training_y = "Normalized_" + filename_training_y


np.save(filename_test_x, TEST_X)
np.save(filename_test_y, TEST_Y)
np.save(filename_training_x, TRAINING_X)
np.save(filename_training_y, TRAINING_Y)
