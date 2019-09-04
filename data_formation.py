import numpy as np
import os
import system_simulator as ss
import random

lstdr = os.listdir()

# Network_files handling
networkQuery = []
for el in lstdr:
    if "nutrient_plant" in el and "txt" in el and "result" not in el:
        networkQuery.append(el)
networkQuery.sort()

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
        resultFiles.append(el)
resultFiles.sort()

numElements = 0

DATA_tuple = []
for i, el in enumerate(resultFiles):
    txt = open(el)
    print("reading + " + el + "...")
    first_line = txt.readline()
    if not numElements:
        numElements = len(first_line.split(","))

    for line in txt:
        splt = line.strip().split(", ")
        if len(splt) != numElements:
            print(el + " has wrong line")
            print(line)
            continue
        X = INITIAL_STATE - np.asarray(splt, dtype=np.float64)
        Y = byunghyun_coeficients[i]
        DATA_tuple.append((X, Y))

print("Total " + str(len(DATA_tuple) + " date prepared."))


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
print("The size of Test Data is : " + str(len(Test_X)))

np.save("test_X.npy", TEST_X)
np.save("test_Y.npy", TEST_Y)
np.save("training_X.npy", TRAINING_X)
np.save("training_Y.npy", TRAINING_Y)

print("DATA SAVED")