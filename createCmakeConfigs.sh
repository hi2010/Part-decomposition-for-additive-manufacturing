#!/bin/bash
cmake -S . -B build/debug -DCMAKE_BUILD_TYPE=DEBUG
cmake -S . -B build/release -DCMAKE_BUILD_TYPE=RELEASE
