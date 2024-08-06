#!/bin/bash

cmake -S./ -B./build -G Ninja -DCMAKE_C_COMPILER=gcc-14
cmake --build ./build
