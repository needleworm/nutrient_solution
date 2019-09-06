import numpy as np
import os
import system_simulator as ss
import random

resultdir = "nutrient_plant_interaction_results/"
lstdr = os.listdir(resultdir)

# Network_files handling
networkQuery = []
for el in lstdr:
    if "nutrient_plant" in el and "txt" in el and "result" not in el:
        networkQuery.append(resultdir + el)
networkQuery.sort()

print(str(len(networkQuery)) + " network query files found.")

byunghyun_coeficients = []
for el in networkQuery:
    temp_network = ss.Network(el)
    byunghyun_coeficients.append(temp_network.byunghyun_coefficients())

INITIAL_STATE = temp_network.Xs
del(temp_network)

# Result_files handling
resultFiles = []
for el in lstdr:
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
        X = INITIAL_STATE - np.asarray(splt, dtype=np.float64)
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

np.save("test_X.npy", TEST_X)
np.save("test_Y.npy", TEST_Y)
np.save("training_X.npy", TRAINING_X)
np.save("training_Y.npy", TRAINING_Y)

print("DATA SAVED")
