# ArCubeOptSim
Geant4 code for ArgonCube optical simulations


# Instructions for installation on lheppc20 machine:
1. Clone this repository: git clone https://github.com/LHEP-neutrino/ArCubeOptSim.git
2. install cmake (version 4.0.3 worked for me)
3. add the following lines to the bashrc:
    export PATH="`pwd`/cmake-4.0.3-linux-x86_64/bin:$PATH"
    source /home/livio/software/root_v6.10.08/bin/thisroot.sh
    source /home/livio/software/geant4/bin/geant4.sh
# Only if you need to edit the source code, or the binary contained in the repository does not work:
4. Open cmake (type cmake-gui in the terminal)
5. Drag-and-drop the CmakeLists.txt from the ArCubeOptSim folder into the cmake window
6. Put the following options in cmake: 
    Geant4_DIR  /home/livio/software/geant4/lib/Geant4-10.4.3
    ROOT_DIR    /home/livio/software/root_v6.10.08/cmake
7. Specify the bin/ folder within ArCubeOptSim as the location to build the binaries. The compiler should be found automatically. If it is not, you need to specify in the cmake options (tick the 'Advanced' box to see it):
    CMAKE_CXX_COMPILER  /usr/bin/c++
    CMAKE_C_COMPILER    /usr/bin/cc
8. In cmake: click 'Configure', then 'Generate', then go into the terminal, navigate to ArCubeOptSim/bin and type 'make -j6'. If this all finishes successfully, you are ready to edit the source code and compile your changes.

# HOW TO RUN:
An example of a command to run the simulation is contained in 'runPDEsim.sh'. Note that you have to change './build/ArgonCubeOptPh' to './bin/ArgonCubeOptPh'. -g specifies the geometry, -p and -m the macros, and -o the name of the output file. There is a preinit macro that controls things like verbosity, step length and TPB layer thickness, and a main macro that controls properties of the generated particles, loads the optical properties of the materials (/argoncube/detector/optical/loadOptSett) and defines the volumes that are sensitive to photon hits (/argoncube/analysis/DefOptSD, valid options are 'volSiPM_Sens_PV' for the ArcLights and 'volSiPM_LCM_PV' for the LCMs). At the end, '/run/beamOn' specifies the amount of runs (the amount of particles per run is specified above at '/argoncube/gun/primaryNb').

# How to add new materials:
This is my best guess, no guarantee for this to work! 
For example, PEN instead of TPB: First, add a new entry in the json file you want to use in the simulation (e.g. resources/OptSim_acl.json for the ArcLight) that has the same properties as the TPB entry (ABSLENGTH, WLSABSLENGTH etc.), then create corresponding .dat files and fill them with numbers from the literature.

# How to edit the geometries:
This is my best guess, no guarantee for this to work! 
For example, replace the LCMs in the full module with ArcLights: in the corresponding gdml file (resources/gdml/module0.gdml), you can find the definitions and locations of all the individual components of the module. In this case, you could simply replace the volumeref "VolLCMPlane" with "volArCLight" wherever it occurs.