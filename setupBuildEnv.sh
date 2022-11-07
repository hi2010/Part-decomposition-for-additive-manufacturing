#!/bin/bash
./installDependencies.sh
./rebuildlibs.sh
./createCmakeConfigs.sh
./buildDebug.sh
./buildRelease.sh
