import numpy as np
from subprocess import run
import os
import sys

# modify stability.txt
def modifFile(dt, dx, end_time, Du, Dv, f, k, r_treshold, sampling, param_file):
    nbParam = 9
    param = [0] * nbParam
    param[0] = dt
    param[1] = dx
    param[2] = end_time
    param[3] = Du
    param[4] = Dv
    param[5] = f
    param[6] = k
    param[7] = r_treshold
    param[8] = int(sampling)
    with open(param_file, "w") as file:
        for i in range(nbParam):
            file.write(str(param[i]) + '\n')

dt = 1e-2
dx_240 = 1e-5 # dx for 10 nodes and 24 threads/process, total time of ~3-4 min
end_time = 1
Du = 1e-5
Dv = 4e-6
f = 0.04
k = 0.06
r_treshold = 1e-5
sampling = 0

nProcess = int(sys.argv[1])
param_file = sys.argv[2]
# nThreads = np.array([24, 20, 16, 12, 8, 4, 2, 1])
nThreads = np.array([4, 2, 1])

dx = dx_240 / np.sqrt(nThreads * nProcess / 240)

for i in range(len(nThreads)):
        modifFile(dt, dx[i], end_time, Du, Dv, f, k, r_treshold, sampling, param_file)
        run(["echo", "*** running with " + str(nThreads[i]) + " threads/process, " + str(nProcess) + " processes and dx = "  + str(dx[i]) + " ***"])
        os.environ["OMP_NUM_THREADS"] = str(nThreads[i])
        run("./run.sh")
print("Test finished wih success")
