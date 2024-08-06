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

# Set a flag for whether to build a cpython shared library
option(BUILD_CPYTHON_SHARED_LIB "Build a shared library for use with CPython" ON)

set(Fortran_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Set source files
set(MAIN_SOURCE
  libbasindiffusion.f90
)
set(AUX_SOURCE
  linalg_tools.f90
  linalg_debug.f90
  oceanography_tools.f90
  constants.f90
)

# Create a list with the resolved paths to the auxiliary source files
set(AUX_SOURCE_FULL "")

foreach(source_file ${AUX_SOURCE})
  list(APPEND AUX_SOURCE_FULL ${Fortran_SOURCE_DIR}/${source_file})
endforeach()

# Set temporary directory for preprocessed files
set(CMAKE_PREPROCESSED_DIR ${CMAKE_BINARY_DIR}/preprocessed)

# Append "_preprocessed" to the main source file name
string(REGEX REPLACE "\\.f90$" "_preprocessed.f90" MAIN_SOURCE_PREPROCESSED ${MAIN_SOURCE})  

# Create a target which produces the source *.f90 files after the preprocessor output
add_custom_target(
  main_preprocessed
  DEPENDS ${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}
)

# Print out the preprocessed source files
foreach(source_file ${SOURCE_PREPROCESSED})
  message(STATUS "Preprocessed source file: ${source_file}")
endforeach()

# Add a custom command to run the preprocessor
add_custom_command(
  OUTPUT ${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}
  COMMAND ${CMAKE_Fortran_COMPILER}
    ${CMAKE_Fortran_FLAGS}
    -E
    -cpp
    -P
    ${Fortran_SOURCE_DIR}/${MAIN_SOURCE}
    > ${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}
  DEPENDS "${Fortran_SOURCE_DIR}/${MAIN_SOURCE}"
)

message(STATUS ${PREPROCESS_CMD})

include(CMakePrintHelpers)

if(BUILD_CPYTHON_SHARED_LIB)
  # Grab Python, 3.11 or newer
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
  cmake_print_variables(Python_INCLUDE_DIRS)
  cmake_print_variables(F2PY_INCLUDE_DIR)
  cmake_print_variables(Python_NumPy_INCLUDE_DIRS)

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
                     "${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}"
                     --lower # Important
                    #  --verbose
    DEPENDS ${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED} # Fortran source
  )

  # Set up target
  Python_add_library(${CMAKE_PROJECT_NAME} MODULE WITH_SOABI
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" # Generated
    "${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
    "${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}" # Fortran source(s)
    "${CMAKE_CURRENT_BINARY_DIR}/${f2py_fortran_wrapper}" # Fortran wrapper
    ${AUX_SOURCE_FULL} # Aux source
  )

  set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "")

  # Depend on sources
  target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE Python::NumPy)
  add_dependencies(${CMAKE_PROJECT_NAME} genpyf)
  target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE "${F2PY_INCLUDE_DIR}")
else()
  add_library(${CMAKE_PROJECT_NAME} SHARED
    ${Fortran_SOURCE_DIR}/${MAIN_SOURCE}
    ${Fortran_SOURCE_DIR}/${AUX_SOURCE}
  )
endif()

set(MODULE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/mod)
set(CMAKE_Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})
set_directory_properties(PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})

set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIRECTORY})

## GSW
# Add the directory for GSW
add_subdirectory(extern/gsw)

## LAPACK
find_package(LAPACK REQUIRED)
include_directories(${LAPACK_INCLUDE_DIRS})

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES} gsw_static)

