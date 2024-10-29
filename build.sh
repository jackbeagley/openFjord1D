#!/bin/bash

cmake -S./src -B./build -DGSW_DIR=../extern/gsw
cmake --build ./build
