import os
import shutil
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
n = np.array(vox_tab.readline().split()).astype(np.int)

print '\n=============================='
print 'number of voxels: %d' % np.prod(nvox)
print 'number of events: %d' % n[0]
print 'number of photons/evt: %d' % n[1]
print '==============================\n'

for vox in range(np.prod(nvox)):

    if (vox == np.prod(nvox)):
        break

    #read voxel information (no. & center coordinates)
    line = np.array(vox_tab.readline().split()).astype(np.float)
    nv = line[0] # voxel number
    cv = line[1:4] # voxel coordinates origin

    #only write macros asked for
    if nv not in range(vox1,vox2):
        continue

    if (vox%1 == 0):
        print "  writing macro no. %d..." % vox

    #create new macro file
    if usro:
        shutil.copy2('macros/OptSim_template_usr.mac','/output/macro_files/OptSim_%08d.mac' % nv)
    else:
        shutil.copy2('macros/OptSim_template.mac','/output/macro_files/OptSim_%08d.mac' % nv)

    macro = open('/output/macro_files/OptSim_%08d.mac' % nv, 'a')

    #write GEANT4 commands into macro file
    macro.write("/argoncube/gun/center %f %f %f mm\n" % (cv[0],cv[1],cv[2]))
    macro.write("/argoncube/gun/halfx " + str(vdim[0]/2.) + " mm\n")
    macro.write("/argoncube/gun/halfy " + str(vdim[1]/2.) + " mm\n")
    macro.write("/argoncube/gun/halfz " + str(vdim[2]/2.) + " mm\n")

    #number of photons per event
    macro.write("/argoncube/gun/primaryNb %d\n" % n[1])

    #number of events
    macro.write("/run/beamOn %d" % n[0])

    macro.close()
