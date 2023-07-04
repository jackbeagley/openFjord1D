cmake_minimum_required(VERSION 3.18) # Needed to avoid requiring embedded Python libs too

project(basin_diffusion
  VERSION 1.0
  DESCRIPTION "Fjord Basin Diffusion tools"
  LANGUAGES C Fortran
)

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
message(
  FATAL_ERROR
  "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Set source files
set(SOURCE
  src/libbasindiffusion.f90
  src/debug_tools.f90
)

# Grab Python, 3.8 or newer
find_package(Python 3.11 REQUIRED
  COMPONENTS Interpreter Development.Module NumPy)

# Grab the variables from a local Python installation
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Print out the discovered paths
include(CMakePrintHelpers)
cmake_print_variables(Python_INCLUDE_DIRS)
cmake_print_variables(F2PY_INCLUDE_DIR)
cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

# Set up the library
# add_library(${CMAKE_PROJECT_NAME} SHARED
#   src/libbasindiffusion.f90
# )

# Common variables
set(f2py_module_name "${CMAKE_PROJECT_NAME}")
set(fortran_src_file "${CMAKE_SOURCE_DIR}/src/libbasindiffusion.f90")
set(f2py_module_c "${f2py_module_name}module.c")
set(f2py_fortran_wrapper "${f2py_module_name}-f2pywrappers2.f90")

# Generate sources
add_custom_target(
  genpyf
  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}"
)
add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" "${CMAKE_CURRENT_BINARY_DIR}/${f2py_fortran_wrapper}"
  COMMAND ${Python_EXECUTABLE}  -m "numpy.f2py"
                   -m "${f2py_module_name}"
                   "${fortran_src_file}"
                   --lower # Important
                  #  --verbose
  DEPENDS ${fortran_src_file} # Fortran source
)

# Set up target
Python_add_library(${CMAKE_PROJECT_NAME} MODULE WITH_SOABI
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" # Generated
  "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
  "${SOURCE}" # Fortran source(s)
  "${CMAKE_CURRENT_BINARY_DIR}/${f2py_fortran_wrapper}" # Fortran wrapper
)

set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "")

# Depend on sources
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE Python::NumPy)
add_dependencies(${CMAKE_PROJECT_NAME} genpyf)
target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE "${F2PY_INCLUDE_DIR}")

set(MODULE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})
set_directory_properties(PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})

set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})

# Add the directory for GSW
add_subdirectory(extern/gsw)

get_property(subdirectory_targets DIRECTORY extern/gsw PROPERTY BUILDSYSTEM_TARGETS)

foreach(target ${subdirectory_targets})
   set_target_properties(${target} PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})
   message(STATUS "Setting Fortran module directory for target ${target}: ${MODULE_OUTPUT_DIRECTORY}")
endforeach()

# Link the shared library to the static library
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE gsw_static)

find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})
target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES})

