#! /bin/bash

cd ./external/vtk
rm -r ./build
mkdir ./build
cd ./build
cmake ..
cmake --build . -j 20

cd ../../VHACD_Performance/VHACD_Lib
rm -r ./build
mkdir ./build
cd ./build
cmake ..
cmake --build . -j 20
