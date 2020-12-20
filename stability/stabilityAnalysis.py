import numpy as np
import os


path = './results/'
filesNameArray = os.listdir(path)

nbFiles = len(filesNameArray)
isStable = np.ndarray(shape=(nbFiles,), dtype='bool')
N = np.ndarray(shape=(nbFiles,), dtype='int')
CourantNumber = np.ndarray(shape=(nbFiles,), dtype='double')

i = 0
for i in range(nbFiles):
    fileName = filesNameArray[i]
    tmp = fileName.split('_')
    N[i] = tmp[1]
    CourantNumber[i] = tmp[2]
    with open(path + fileName, "r") as file:
        n = np.fromfile(file, dtype=np.int32, count=1)
        arr = np.fromfile(file, dtype=np.double, count=n[0]**2, offset=4)
        print(arr)
        if np.isnan(arr).any():
            isStable[i] = False
        else:
            isStable[i] = True

print(isStable)
print(N)
print(CourantNumber)
