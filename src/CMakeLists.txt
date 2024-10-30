cmake_minimum_required(VERSION 3.18...3.22)# Needed to avoid requiring embedded Python libs too

project(basin_diffusion
  VERSION 1.0
  DESCRIPTION "Fjord Basin Diffusion tools"
  LANGUAGES C Fortran
)

if (NOT SKBUILD)
	message(STATUS "Build not controlled by scikit-build")
else()
	message(STATUS "Build controlled by scikit-build")
endif()

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
message(
  FATAL_ERROR
  "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Set a flag for whether to build a cpython shared library
option(BUILD_CPYTHON_SHARED_LIB "Build a shared library for use with CPython" ON)

set(Fortran_SOURCE_DIR ${CMAKE_SOURCE_DIR})

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
file(MAKE_DIRECTORY ${CMAKE_PREPROCESSED_DIR})

# Append "_preprocessed" to the main source file name
string(REGEX REPLACE "\\.f90$" "_preprocessed.f90" MAIN_SOURCE_PREPROCESSED ${MAIN_SOURCE})  

# Create a target which produces the source *.f90 files after the preprocessor output
add_custom_target(
  main_preprocessed
  DEPENDS ${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}
)

# Print out the preprocessed source files
#foreach(source_file ${SOURCE_PREPROCESSED})
#  message(STATUS "Preprocessed source file: ${source_file}")
#endforeach()

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
	# Ensure scikit-build modules
	if (NOT SKBUILD)
  		find_package(Python 3.8 REQUIRED)
	
		# Set the Python executable for compatibility
		set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})

  		# Kanged --> https://github.com/Kitware/torch_liberator/blob/master/CMakeLists.txt
		# If skbuild is not the driver; include its utilities in CMAKE_MODULE_PATH
		execute_process(
		    COMMAND "${PYTHON_EXECUTABLE}"
		    -c "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
		    OUTPUT_VARIABLE SKBLD_DIR
		    OUTPUT_STRIP_TRAILING_WHITESPACE
  		)
  		list(APPEND CMAKE_MODULE_PATH "${SKBLD_DIR}/resources/cmake")
  		message(STATUS "Looking in ${SKBLD_DIR}/resources/cmake for CMake modules")
	endif()

	# scikit-build style includes
	find_package(PythonExtensions REQUIRED) # for ${PYTHON_EXTENSION_MODULE_SUFFIX}

	# Grab the variables from a local Python installation
	# NumPy headers
	execute_process(
	  COMMAND "${PYTHON_EXECUTABLE}"
	  -c "import numpy; print(numpy.get_include())"
	  OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
	  OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	
	# F2PY headers
	execute_process(
		COMMAND "${PYTHON_EXECUTABLE}"
		-c "import numpy.f2py; print(numpy.f2py.get_include())"
		OUTPUT_VARIABLE F2PY_INCLUDE_DIR
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

  	# Print out the discovered paths
	cmake_print_variables(PYTHON_INCLUDE_DIRS)
  	cmake_print_variables(F2PY_INCLUDE_DIR)
  	cmake_print_variables(NumPy_INCLUDE_DIRS)

  	set(f2py_module_name "${CMAKE_PROJECT_NAME}")
  	set(fortran_src_file "${CMAKE_SOURCE_DIR}/libbasindiffusion.f90")
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
			--quiet
		#--verbose
		DEPENDS ${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED} # Fortran source
	)

	# Set up target
	add_library(${CMAKE_PROJECT_NAME} MODULE
		"${CMAKE_CURRENT_BINARY_DIR}/${f2py_module_c}" # Generated
		"${F2PY_INCLUDE_DIR}/fortranobject.c" # From NumPy
		"${CMAKE_PREPROCESSED_DIR}/${MAIN_SOURCE_PREPROCESSED}" # Fortran source(s)
		"${CMAKE_CURRENT_BINARY_DIR}/${f2py_fortran_wrapper}" # Fortran wrapper
		${AUX_SOURCE_FULL} # Aux source
	)

	target_include_directories(${CMAKE_PROJECT_NAME} PUBLIC
		${F2PY_INCLUDE_DIR}
		${NumPy_INCLUDE_DIRS}
		${PYTHON_INCLUDE_DIRS})
	set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES SUFFIX "${PYTHON_EXTENSION_MODULE_SUFFIX}")
	set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES PREFIX "")
	set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "")

	# Linker fixes
	if (UNIX)
	  if (APPLE)
	    set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES
	    LINK_FLAGS  '-Wl,-dylib,-undefined,dynamic_lookup')
	  else()
	    set_target_properties(${CMAKE_PROJECT_NAME} PROPERTIES
	  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
	  endif()
	endif()

	# Depend on sources
	add_dependencies(${CMAKE_PROJECT_NAME} genpyf)
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
add_subdirectory(${EXTERN_DIR}/gsw gsw)

## LAPACK
find_package(LAPACK)

if (LAPACK_FOUND)
	include_directories(${LAPACK_INCLUDE_DIRS})
else()
	add_subdirectory(${EXTERN_DIR}/lapack lapack)
endif()

target_link_libraries(${CMAKE_PROJECT_NAME} PRIVATE ${LAPACK_LIBRARIES} gsw_static)

