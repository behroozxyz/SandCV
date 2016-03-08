#==============================================================================
# File         : CMakeLists.txt 
# Copyright    : (C)opyright 2013
# Authors      : Behrooz Hashemian
# Last modified: Fri Jan 30, 2015  10:46AM
#============================================================================*/
cmake_minimum_required(VERSION 2.8)

#------------------------------------------------------------------------------
# SandCV version number.  An even minor number corresponds to releases.
#------------------------------------------------------------------------------
#set(SandCV_VERSION_MAJOR "0")
#set(SandCV_VERSION_MINOR "1")
#set(SandCV_VERSION_PATCH "0")
#set(SandCV_VERSION_STRING "${SandCV_VERSION_MAJOR}.${SandCV_VERSION_MINOR}.${SandCV_VERSION_PATCH}")

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------
option(BUILD_DOCUMENTATION "Build Doxygen documentation"         ON )
option(BUILD_TESTING       "Build SandCV test applications"      OFF)
option(BUILD_EXAMPLES      "Build SandCV Examples"               ON )
option(BUILD_SHARED_LIBS   "Build SandCV with shared libraries"  ON )
option(USE_MPI             "Use Message Passing Interface (MPI)" OFF)
option(USE_MKL             "Use Intel Math Kernel Library"       ON )
option(USE_OPENMP          "Use compiler C/C++ with OpenMP"      ON )
option(USE_MATLAB          "Use Matlab data format (.mat)"       ON )
option(USE_SWIG            "Use SWIG to build TCL interface"     OFF)
option(USE_INTEL           "Use Intel Compiler C/C++"            ON )
option(USE_Boost           "Use Bosst C++ Library"               OFF)

#------------------------------------------------------------------------------
# Settings
#------------------------------------------------------------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build Type")
set(CMAKE_INSTALL_PREFIX "$ENV{HOME}/local" CACHE STRING "Install directory prefix")

# Define settings for Intel Compiler
if(USE_INTEL)
    find_path(USE_INTEL_ROOT icpc)
    set(CMAKE_CXX_COMPILER "${USE_INTEL_ROOT}/icpc" CACHE STRING "Using Intel C++ Compiler")
    mark_as_advanced(USE_INTEL_ROOT)
else(USE_INTEL)
    find_path(GNU_ROOT g++)
    set(CMAKE_CXX_COMPILER "${GNU_ROOT}/g++"  CACHE STRING "Using GNU C++ Compiler")
    mark_as_advanced(GNU_ROOT)
endif(USE_INTEL)

# Define settings for Intel or GCC Compiler in order to use OpenMP
if(USE_OPENMP)
    if(USE_INTEL)
        set(CMAKE_CXX_FLAGS "-openmp" CACHE STRING "Flag needed to activate OpenMP pragmas recognition")
    else(USE_INTEL)
        if(UNIX AND CMAKE_COMPILER_IS_GNUCXX)
            set(CMAKE_CXX_FLAGS "-fopenmp" CACHE STRING "Flag needed to activate OpenMP pragmas recognition")
        endif(UNIX AND CMAKE_COMPILER_IS_GNUCXX)
    endif(USE_INTEL)
endif(USE_OPENMP)

# C++ Compiler options
set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g"   CACHE STRING "Flag to generate debug builds"  )
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG" CACHE STRING "Flag to generate release builds")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG" CACHE STRING "Flag to generate minsize release builds")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g"       CACHE STRING "Flag to generate release builds with debug information")

#------------------------------------------------------------------------------
# The cache variables have to be set before project(),
# to have an overridable predefined value without FORCE!
#------------------------------------------------------------------------------
# PROJECT
#------------------------------------------------------------------------------
project(SandCV CXX)
set(CMAKE_MACOSX_RPATH 1)

#------------------------------------------------------------------------------
# Find required Libraries
#------------------------------------------------------------------------------
# Intel MKL
if(USE_MKL)
    find_package(MKL REQUIRED)
    include_directories(${MKL_INCLUDE_DIRS})
    if (NOT USE_INTEL AND USE_OPENMP)
        set(MKL_LIBRARIES ${MKL_LIBRARIES} mkl_gnu_thread )
    endif (NOT USE_INTEL AND USE_OPENMP)
endif(USE_MKL)

# Boost C++ Library
if(USE_Boost)
    find_package(Boost REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
endif(USE_Boost)


# MATLAB
if(USE_MATLAB)
    find_package(MATLAB REQUIRED)
    include_directories(${MATLAB_INCLUDE_DIRS})
endif(USE_MATLAB)

# Swig
if(USE_SWIG)
    find_package(SWIG REQUIRED)
    include(${SWIG_USE_FILE})
    #  find_package(TCL REQUIRED)
    #  include_directories(${TCL_INCLUDE_PATH})
    find_package(TCL)
    if(TCL_FOUND)
        include_directories(${TCL_INCLUDE_PATH})
    else(TCL_FOUND)
        message(FATAL_ERROR "Could not find TCL.")
    endif(TCL_FOUND)
endif(USE_SWIG)

# Zlib
#set(ZLIB_LIBRARIES z)

#------------------------------------------------------------------------------
# Compilation settings
#------------------------------------------------------------------------------
add_definitions(-std=c++11) # Enable C++11 features in the compiler
if(APPLE)
    add_definitions(-stdlib=libc++) # Enable C++11 features in the compiler
endif(APPLE)
add_definitions(-Wall)      # All warnings
#add_definitions(-ansi)     # Tells the compiler to implement the ANSI language option
if(USE_OPENMP)
    add_definitions(-openmp_report=1) # {0|1|2} Control the OpenMP parallelizer diagnostic level
endif(USE_OPENMP)

#------------------------------------------------------------------------------
# Output directories
#------------------------------------------------------------------------------
set(LIBRARY_OUTPUT_PATH    ${SandCV_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
set(EXECUTABLE_OUTPUT_PATH ${SandCV_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

#link_directories(${SandCV_BINARY_DIR})

#------------------------------------------------------------------------------
# SandCV Configuration files 
#------------------------------------------------------------------------------
configure_file( ${SandCV_SOURCE_DIR}/SandCVConfig.h.in
        ${SandCV_BINARY_DIR}/SandCVConfig.h   )
# Include the directory of SandCVConfig.h
include_directories(
        ${SandCV_BINARY_DIR}
        ${SandCV_SOURCE_DIR}/Modules
        ${SandCV_SOURCE_DIR}/Modules/Alignment
        ${SandCV_SOURCE_DIR}/Modules/ColVars
        ${SandCV_SOURCE_DIR}/Modules/IO
        ${SandCV_SOURCE_DIR}/Modules/Parametrization
        ${SandCV_SOURCE_DIR}/Modules/Projection
        ${SandCV_SOURCE_DIR}/Modules/Searcher
        ${SandCV_SOURCE_DIR}/Modules/Searcher/ANN
        ${SandCV_SOURCE_DIR}/Modules/Simulation
)

#------------------------------------------------------------------------------
# Add subdirectories
#------------------------------------------------------------------------------
add_subdirectory(Modules)
if(USE_SWIG)
    include_directories( ${SandCV_SOURCE_DIR}/Interfaces )
    add_subdirectory(Interfaces)
endif(USE_SWIG)

if(BUILD_TESTING)
    add_subdirectory(Testing)
endif(BUILD_TESTING)

if(BUILD_EXAMPLES)
    add_subdirectory(Examples)
endif(BUILD_EXAMPLES)

#------------------------------------------------------------------------------
# Documentation configure files
#------------------------------------------------------------------------------
if(BUILD_DOCUMENTATION)
    # Doxygen
    find_package(Doxygen REQUIRED)
    configure_file(${SandCV_SOURCE_DIR}/Documentation/Doxygen/Doxyfile.in
            ${SandCV_BINARY_DIR}/Documentation/Doxyfile)
endif(BUILD_DOCUMENTATION)

#------------------------------------------------------------------------------
# SandCV installation structure
#------------------------------------------------------------------------------
if(NOT SandCV_INSTALL_RUNTIME_DIR)
    set(SandCV_INSTALL_RUNTIME_DIR bin)
endif()
if(NOT SandCV_INSTALL_LIBRARY_DIR)
    set(SandCV_INSTALL_LIBRARY_DIR lib/SandCV)
endif()
if(NOT SandCV_INSTALL_INCLUDE_DIR)
    set(SandCV_INSTALL_INCLUDE_DIR include/SandCV-${SandCV_MAJOR_VERSION}.${SandCV_MINOR_VERSION})
endif()
if(NOT SandCV_INSTALL_DATA_DIR)
    set(SandCV_INSTALL_DATA_DIR share/SandCV-${SandCV_MAJOR_VERSION}.${SandCV_MINOR_VERSION})
endif()
if(NOT SandCV_INSTALL_DOC_DIR)
    set(SandCV_INSTALL_DOC_DIR share/doc/SandCV-${SandCV_MAJOR_VERSION}.${SandCV_MINOR_VERSION})
endif()
if(NOT SandCV_INSTALL_DOXYGEN_DIR)
    set(SandCV_INSTALL_DOXYGEN_DIR ${SandCV_INSTALL_DOC_DIR}/doxygen)
endif()

#-----------------------------------------------------------------------------
# Mac platform considerations
#-----------------------------------------------------------------------------
if(APPLE)
    mark_as_advanced(
            CMAKE_OSX_ARCHITECTURES
            CMAKE_OSX_DEPLOYMENT_TARGET
            CMAKE_OSX_SYSROOT
    )
endif(APPLE)

