#-----------------------------------------------------------------------------
# FindMKL.cmake
#
# Module looks for MKL libraries and includes' files.
#
# Defines:
#  MKL_INCLUDE_DIRS : include path
#  MKL_LIBRARIES    : required libraries
#  MKL_RT_LIBRARY   : path to libmkl_rt.so which is single dynamic library and
#                     is supported since MKL 10.3
#
#-----------------------------------------------------------------------------

# The MKL root directory.
if(DEFINED ENV{MKLROOT})
  set(MKL_ROOT  $ENV{MKLROOT})
else(DEFINED ENV{MKLROOT})
  set(MKL_ROOT  /opt/intel/mkl)
  message(WARNING "MKLROOT environment vriable is not set. ${MKL_ROOT} is used.")
endif(DEFINED ENV{MKLROOT})

# The MKL include file directories.
set(MKL_INCLUDE_DIRS ${MKL_ROOT}/include)

# The MKL lib library directories.
if(UNIX)
  set(MKL_LIBRARY_DIRS ${MKL_ROOT}/lib/intel64)
endif(UNIX)
if(APPLE)
  set(MKL_LIBRARY_DIRS ${MKL_ROOT}/lib)
endif(APPLE)

# This find the required libraries
find_library(MKL_RT_LIBRARY         mkl_rt         PATHS ${MKL_LIBRARY_DIRS} NO_SYSTEM_ENVIRONMENT_PATH)
#find_library(MKL_CORE_LIBRARY         mkl_core         HINTS ${MKL_LIBRARY_DIRS})
#find_library(MKL_INTEL_LP64_LIBRARY   mkl_intel_lp64   HINTS ${MKL_LIBRARY_DIRS})
#find_library(MKL_INTEL_THREAD_LIBRARY mkl_intel_thread HINTS ${MKL_LIBRARY_DIRS})

set(MKL_LIBRARIES
  ${MKL_RT_LIBRARY}
#  ${MKL_CORE_LIBRARY}
#  ${MKL_INTEL_LP64_LIBRARY}
#  ${MKL_INTEL_THREAD_LIBRARY}
  )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL DEFAULT_MSG
  MKL_INCLUDE_DIRS 
  MKL_LIBRARIES)

mark_as_advanced(
  MKL_ROOT
  MKL_RT_LIBRARY
#  MKL_CORE_LIBRARY
#  MKL_INTEL_LP64_LIBRARY
#  MKL_INTEL_THREAD_LIBRARY
  )
