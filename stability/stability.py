import numpy as np
from subprocess import run
from os import rename


# modify stability.txt
def modifFile(dt, dx, end_time, Du, Dv, f, k, r_treshold, sampling):
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
    with open("stability.txt", "w") as file:
        for i in range(nbParam):
            file.write(str(param[i]) + '\n')

def retrieveFile(N, CourantNumber, sampling):
    rename(r'./results/u_' + str(int(sampling)) + '.dat', r'./results/u_' + str(N) + '_' + str(CourantNumber) + '_.dat')
    rename(r'./results/v_' + str(int(sampling)) + '.dat', r'./results/v_' + str(N) + '_' + str(CourantNumber) + '_.dat')

end_time = 1e5
Du = 1e-5
Dv = 2e-6
f = 0.04
k = 0.06
r_treshold = 1e-5
sampling = 0

N = np.array([400]) # , 600, 800, 1000]) #, 100, 75, 50, 25])
# N = np.array([100])
# CourantNumber = np.array([0.1])
CourantNumber = np.array([0.2])
dx = 1/N

dt = np.ndarray(shape=(len(dx), len(CourantNumber)))
sampling = np.ndarray(shape=(len(dx), len(CourantNumber)))
for i in range(len(dx)):
    dt[i, :] = CourantNumber * dx[i]**2 / np.max([Du, Dv])

sampling = end_time / dt
for i in range(len(dx)):
    for j in range(len(CourantNumber)):
        # run(["echo", "$SLURM_NTASKS"])
        modifFile(dt[i, j], dx[i], end_time, Du, Dv, f, k, r_treshold, sampling[i, j])
        run(["echo", "*** running with CourantNumber = " + str(CourantNumber[j]) + ", and N = " + str(N[i]) + ", s = " + str(sampling[i, j]) + " ***"])
        run("./run.sh")
        retrieveFile(N[i], CourantNumber[j], sampling[i, j])
print("files moved with success")
