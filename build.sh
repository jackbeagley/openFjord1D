#!/bin/bash

cmake -S./src -B./build -G Ninja -DEXTERN_DIR=../extern
cmake --build ./build
