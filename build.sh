#!/bin/bash

cmake -S./src -B./build -G Ninja -DGSW_PATH=../extern/gsw
cmake --build ./build
