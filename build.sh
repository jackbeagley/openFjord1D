#!/bin/bash

cmake -S./openFjord1D -B./build -G Ninja -DEXTERN_DIR=../extern
cmake --build ./build
