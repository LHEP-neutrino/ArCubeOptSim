#!/bin/bash


#Uncomment to bind mount ArCubeOptSim project directory for building it
docker run -it --rm --name geant4_container_build -v $(pwd)/build:/tmp/build -v $(pwd)/repo:/tmp/code smflment/geant4:10.5.1 bash



