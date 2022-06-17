import os
import subprocess
import numpy as np

#read environent variables
usrg = int(os.environ['USRG'])

#read voxel table
if usrg:
    vox_tab = open('/input/OptSim_LUT_voxel_table.txt', 'r')
else:
    vox_tab = open('OptSim_LUT_voxel_table.txt', 'r')

cmin = np.array(vox_tab.readline().split()).astype(np.float)
cmax = np.array(vox_tab.readline().split()).astype(np.float)
vdim = np.array(vox_tab.readline().split()).astype(np.float)
nvox = np.array(vox_tab.readline().split()).astype(np.int)
print '\n------------------------------'

vox_per_bunch = 100

#only LUT sim files as bunches of 100
for bunch in range(np.prod(nvox)/vox_per_bunch+1):

    vox_min = bunch*vox_per_bunch
    if (bunch==np.prod(nvox)/vox_per_bunch):
        vox_max = np.prod(nvox)%vox_min
    else:
        vox_max = vox_per_bunch

    print "\n===> processing voxels %d to %d..." % (vox_min, vox_min+vox_max-1)

    root_call = "OptSim_createLUT.C(%d)" % bunch

    subprocess.call(["/Software/root_v6.22.02/bin/root.exe","-l","-q",root_call])

subprocess.call(["rm", "/output/OptSim_LUT_ArgonCube2x2.root"])
subprocess.call('/Software/root_v6.22.02/bin/hadd /output/OptSim_LUT_ArgonCube2x2.root /output/*.root', shell=True)
