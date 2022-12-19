import os
import subprocess
import numpy as np

#read environent variables
usrg = int(os.environ['USRG'])
sipm_eff = float(os.environ['SIPMEFF'])
start_vox = int(os.environ['STARTVOX'])
stop_vox = int(os.environ['STOPVOX'])

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
if (stop_vox == -1):                #run over all voxels
    stop_vox = np.prod(nvox)

#only LUT sim files as bunches of 100
for bunch in range(np.prod(nvox)/vox_per_bunch+1):
    vox_min = bunch*vox_per_bunch
    if (vox_min<start_vox):
	continue
    if (bunch==np.prod(nvox)/vox_per_bunch):
        vox_max = np.prod(nvox)%vox_min
    else:
        vox_max = vox_per_bunch

    print "\n===> processing voxels %d to %d..." % (vox_min, vox_min+vox_max-1)

    if (vox_min+vox_max >= stop_vox):
        exit(0)

    root_call = "OptSim_createLUT.C(%d,%f)" % (bunch, sipm_eff)

    subprocess.call(["/Software/root_v6.22.02/bin/root.exe","-l","-q",root_call])

subprocess.call(["rm", "/output/OptSim_LUT_ArgonCube2x2.root"])
subprocess.call('/Software/root_v6.22.02/bin/hadd /output/OptSim_LUT_ArgonCube2x2.root /output/*.root', shell=True)
