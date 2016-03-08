#-----------------------------------------------------------------------------
# FindMATLAB.cmake
#
# Module looks for MATLAB libraries and includes' files.
#
# Defines:
#  MATLAB_INCLUDE_DIRS   : include path for mex.h, engine.h
#  MATLAB_LIBRARIES  : required libraries: libmex, etc
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY : path to libmx.lib
#  MATLAB_ENG_LIBRARY: path to libeng.lib
#
#-----------------------------------------------------------------------------

# The MATLAB root directory
if(DEFINED ENV{MATLABROOT})
  set(MATLAB_ROOT  $ENV{MATLABROOT})
else(DEFINED ENV{MATLABROOT})
  if(UNIX)
    set(MATLAB_ROOT /opt/MATLAB/R2013a)
  endif(UNIX)
  if(APPLE)
    set(MATLAB_ROOT /Applications/MATLAB_R2012b.app)
  endif(APPLE)
  message(WARNING "MATLABROOT environment vriable is not set. ${MATLAB_ROOT} is used.")
endif(DEFINED ENV{MATLABROOT})

# The MATLAB lib library directories.
if(UNIX)
  set(MATLAB_LIBRARY_DIRS ${MATLAB_ROOT}/bin/glnxa64)
endif(UNIX)
if(APPLE)
  set(MATLAB_LIBRARY_DIRS ${MATLAB_ROOT}/bin/maci64)
endif(APPLE)

# The MATLAB include file directories.
set(MATLAB_INCLUDE_DIRS ${MATLAB_ROOT}/extern/include)


# This find the required libraries
find_library(MATLAB_MX_LIBRARY  mx  HINTS ${MATLAB_LIBRARY_DIRS})
find_library(MATLAB_MEX_LIBRARY mex HINTS ${MATLAB_LIBRARY_DIRS})
find_library(MATLAB_MAT_LIBRARY mat HINTS ${MATLAB_LIBRARY_DIRS})
find_library(MATLAB_ENG_LIBRARY eng HINTS ${MATLAB_LIBRARY_DIRS})

set(MATLAB_LIBRARIES
  ${MATLAB_MX_LIBRARY}
  ${MATLAB_MEX_LIBRARY}
  ${MATLAB_MAT_LIBRARY}
  ${MATLAB_ENG_LIBRARY}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MATLAB DEFAULT_MSG
  MATLAB_INCLUDE_DIRS
  MATLAB_LIBRARIES)

mark_as_advanced(
  MATLAB_ROOT
  MATLAB_MX_LIBRARY
  MATLAB_MEX_LIBRARY
  MATLAB_MAT_LIBRARY
  MATLAB_ENG_LIBRARY
)
