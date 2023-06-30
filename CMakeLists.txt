cmake_minimum_required(VERSION 3.10)
project(BasinDiffusionLibrary Fortran)

# Set the path to the static library folder
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/libs)

# Set the source files for the Fortran library
set(SOURCES
    src/libbasindiffusion.f90
)

# Add the shared library target
add_library(basin_diffusion SHARED ${SOURCES})

# Set the output name of the shared library
set_target_properties(basin_diffusion PROPERTIES OUTPUT_NAME basin_diffusion)

# Add the directory for GSW
add_subdirectory(extern/gsw)

# Link the shared library to the static library
target_link_libraries(basin_diffusion PRIVATE gsw_static)
