import os
import subprocess
import numpy as np

#read environent variables
usrg = int(os.environ['USRG'])
usro = int(os.environ['USRO'])
vox1 = int(os.environ['VOX1'])
vox2 = int(os.environ['VOX2'])

#read voxel table
if usrg:
    vox_tab = open('/input/OptSim_LUT_voxel_table.txt', 'r')
else:
    vox_tab = open('OptSim_LUT_voxel_table.txt', 'r')

cmin = np.array(vox_tab.readline().split()).astype(np.float)
cmax = np.array(vox_tab.readline().split()).astype(np.float)
vdim = np.array(vox_tab.readline().split()).astype(np.float)
nvox = np.array(vox_tab.readline().split()).astype(np.int)

print '\n------------------------------\n'
#print 'Number of voxels: %d' % np.prod(nvox)

gdml = "../resources/gdml/OptSim.gdml"

if os.path.isfile('/input/preinit.mac'):
    print '\nuser-defined preinit macro (TPC thickness)...'
    preinit = "/input/preinit.mac"
else:
    preinit = "../resources/macros/preinit.mac"

#only process voxels asked for
for vox in range(vox1,vox2):

    if (vox == np.prod(nvox)):
        break

    if (vox%1 == 0):
        print "  running simulation on voxel no. %d..." % vox

    infile  = "/output/macro_files/OptSim_%08d.mac" % vox
    outfile = "/output/root_files/OptSim_%08d.root" % vox
    logfile = "/output/log_files/OptSim_%08d.log" % vox
    errfile = "/output/log_files/OptSim_%08d.txt" % vox

    with open(logfile,"wb") as out, open(errfile, "wb") as err:
        subprocess.call(["../build/ArgonCubeOptPh",
            "-g", gdml,
            "-p", preinit,
            "-m", infile,
            "-o", outfile,
            ], stdout=out, stderr=err)

