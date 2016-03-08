if(USE_SWIG)
  #-----------------------------------------------------------------------------
  # TCL interface for Adaptive Biasing Force with Multi-Pactch SandCV
  #-----------------------------------------------------------------------------
  add_library(MPSandRest MPSandRest.cpp)
  target_link_libraries(MPSandRest SandCV IO)
  #-----------------------------------------------------------------------------
  # TCL interface for Adaptive Biasing Force with Multi-Pactch SandCV
  #-----------------------------------------------------------------------------
  add_library(MPSandABF MPSandABF.cpp)
  target_link_libraries(MPSandABF ABF SandCV IO)
  #-----------------------------------------------------------------------------
  # TCL interface for Adaptive Biasing Force with Puckering Coordinates
  #-----------------------------------------------------------------------------
  add_library(PuckerABF PuckerABF.cpp)
  target_link_libraries(PuckerABF ABF Puckering IO)
  #-----------------------------------------------------------------------------
  # TCL interface for Adaptive Biasing Force with SandCV
  #-----------------------------------------------------------------------------
  add_library(SandABF SandABF.cpp)
  target_link_libraries(SandABF ABF SandCV IO)
  #-----------------------------------------------------------------------------
  # TCL interface for Steered Molecule Dynamics with SandCV
  #-----------------------------------------------------------------------------
  add_library(SandSMD SandSMD.cpp)
  target_link_libraries(SandSMD SMD SandCV IO)

  add_library(SandCMapSMD SandCMapSMD.cpp)
  target_link_libraries(SandCMapSMD SMD SandCVCMap IO)
  #-----------------------------------------------------------------------------
  # TCL interface for Adaptive Biasing Force with Puckering Coordinates
  #-----------------------------------------------------------------------------
  add_library(ZoneLock ZoneLock.cpp)
  target_link_libraries(ZoneLock SandCV IO)

  set(CMAKE_SWIG_FLAGS "-c++")
  set_source_files_properties(Interface.i PROPERTIES CPLUSPLUS ON)
  set_source_files_properties(Interface.i PROPERTIES CMAKE_SWIG_FLAGS "-includeall")
  swig_add_module(Interface tcl Interface.i)
  swig_link_libraries(Interface ${TCL_LIBRARIES} PuckerABF MPSandRest MPSandABF SandABF SandSMD SandCMapSMD ZoneLock)
endif(USE_SWIG)
