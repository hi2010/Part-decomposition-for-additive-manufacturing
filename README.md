# Part decomposition for additive manufacturing, 3d printing

This program uses a genetic algorithm to optimize a decomposition of a 3d model with 3mf format.
If your file is not in the 3mf format, you can convert it to obj format and use obj23mf.
The provided cmd-line interface can be seen by running the program without argument or with -h.

The source code for the decomposition code.
This project uses cmake.
Will install to /bin folder if cmake install is used.

Uses the decomposition process from:
[Part decomposition and evaluation based on standard design guidelines for additive manufacturability and assemblability] (https://doi.org/10.1016/j.addma.2020.101702)
slightly modified.


# build
Open in VisualStudioCode with cmake extension or use cmake --build .. and cmake --install .. from the build folder.
Use start script from /bin Folder to use from commandline (or fix the rpath using the cmake file in ./ folder).

vtk should automatically be downloaded to /external/vtk folder when pulling the repo using git.

To install dependencies on a debian based system (ubuntu, ...) one can use ./installDependencies.sh
To build the dependecie vtk either follow the instructions on https://vtk.org/Wiki/VTK/Configure_and_Build to only build needed parts of the vtk library or run ./buildDependencies.sh to build with wathever defaults the build process uses if no ccmake configuration is done.
If ./buildDependencies.sh fails try to run it again.
Optionally in this file, edit the number of threads (j) for optimum build performance.

To use this from docker f.e. with ubuntu run:
``` bash
docker run -ti --name decomp ubuntu /bin/bash
# or to restart the container after a computer restart or other reason it stopped:
docker start -i decomp
```

To copy folders from the docker container to the host pc (if you need files you build on the container and dont use a mount (f.e. you use windows -> can't use host mount otherwise symlinks will fail).
```bash
# copy decompositionApp folder from decomp container to current dir "."
docker container cp decomp:decompositionApp .
```

cd /path where the project should be downloaded to, f.e.:
``` bash
cd /home
```

install git
``` bash
apt update
apt install git
```

pull the repo (use link from gitlab after creating a ssh key or use https version):
``` bash
git pull git@git.link to the repo
```

Setting up the neccessary dependencies:
``` bash
cd ./ name of the repo
./setupBuildEnv.sh
```
If git submodule in installDependencies fails, download / clone the needed vtk lib files from https://vtk.org/download/ to the /external/vtk folder.
Then you can try to rerun the ./setupBuildEnv.sh or run the needed scripts one by one.
Needed programs should already be installed at this point.

configure and build (for j20 (threads = 20), replace with your thread count (can be retrieved using lscpu) for fastest build times:
(from ./build)
``` bash
cd ./name of the repo
cd ./build
cmake ..
cmake --build . -j 20
```
for release and debug builds createCmakeConfigs.sh (run in setupBuildEnv.sh) should create the folders release and debug in the build folder containing the needed cmake files to build a release or debug version (should use different compiler optimization options if it works).
From tese folder use cmake --build . -j 20 to build the program. Or cmake --install . to install (cp) the files to the /bin folder.
``` bash
cd build/release
cmake --build . -j 20
cmake --install .
cd build/debug
cmake --build . -j 20
```

install to the ./bin folder (from ./build):
``` bash
cmake --install .
```


run on the compiile system:
``` bash
./decomposer add arguments here
```
! This only works if $PWD is /build (if one executes the program from within the /build folder / from a folder where lib3mf.so.2 lies (see ldd))

run from any system:
(in ./bin folder)
``` bash
./runDecomposer.sh add arguments here
```
This is needed becaus the script adds the ./bin/lib folder to the ldlibrary path so the system can find the libs.
This could prbly also be achieved by "fixing" the rpath in cmake to point to the correct folder.
Include paths of dynamic linked libraries can be shown by using ldd ./decomposer.
Using the script works and is flexible, as the lib folder can easily be changed in the script.

## Important notes
If one wants to edit the files on a windows system or start the docker container with the service from a windows system, the lib folder within the /bin install folder needs to get zipped using the procedure described in the readme of the service.
This is important, because otherwise the symlinks of libvtk get destroyed the moment the repo gets pulled to the windows filesystem, at least this happened in one test.

For use in the service, copy the contents of the /bin folder to the appropriate folder in the service's project.

## Setup of the wanted decomposition
Adjust the settings in ./decompositionApp/decomposeConfig.json to your likings
and append another list entry if needed.
Per entry one decomposition process gets started. \
Parameter description:
- costWeights: the weights used for the optimisation of the decomposition. Describes the importance of the quality criteria. For more infos see: https://doi.org/10.1016/j.addma.2020.101702 .
- gaOpts: parameters of the optimization process.
- printerType: irrelevant / only "xyz" is supported but could be used for printer build space description or kinematics.
- printerDims: the target dimensions of the part in mm. The parts units are considered to be in mm.
- configName: the name of this configuration.

If costs for per part budget and decomposition cost calculation should be adjusted: \
Adjust the values in DecompositionService at the beginning or any of the formulas.

## How to compress the build to keep symlinks on linux:
```bash
zip -r --symlinks bin.zip bin
```
With bin beeing the folder containing the code.
Use
```bash
unzip bin.zip
```
to unzip and restore the symlinks.

# used libraries:
* [vtk](https://vtk.org/)
* [lib3mf](https://github.com/3MFConsortium/lib3mf)
* [openGa](https://github.com/Arash-codedev/openGA)
* [v-hacd](https://github.com/kmammou/v-hacd) slightly modified version from performance improvements branch

For licenses see licenses.txt.
Starting the program with --licenses as parameter should also show the licenses.

# build system:
[cmake](https://cmake.org/download/)
