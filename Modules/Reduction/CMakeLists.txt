#------------------------------------------------------------
# Source files for smDimensionReduction
#------------------------------------------------------------
SET (DIMREDUCTION_SOURCES_H     smDimensionReduction.h
                                smEmbeddedSamplePositions.h
                                smClosestPointProjection.h
                                smClosestPointProjectionMD.h
                                smClosestPointProjectionLS.h
                                smPCA.h
                                smTangentProjection.h
                                smCorral.h
                                smContinuousCorral.h
                                smMapping.h
                                smSearching.h
                                smAcceleratedMD.h)

SET (DIMREDUCTION_SOURCES_CPP   smDimensionReduction.cpp
                                smEmbeddedSamplePositions.cpp
                                smClosestPointProjection.cpp
                                smClosestPointProjectionMD.cpp
                                smClosestPointProjectionLS.cpp
                                smPCA.cpp
                                smTangentProjection.cpp
                                smCorral.cpp
                                smContinuousCorral.cpp
                                smMapping.cpp
                                smSearching.cpp
                                smAcceleratedMD.cpp)

IF (SOLMEC_USE_ARPACK AND SOLMEC_USE_MKL)
    SET(DIMREDUCTION_SOURCES_H   ${DIMREDUCTION_SOURCES_H}
                                 smLocalIsometricEmbedding.h
                                 smNonlinearIsometricEmbedding.h
                                 smLLE.h smIsomap.h  )
    SET(DIMREDUCTION_SOURCES_CPP ${DIMREDUCTION_SOURCES_CPP}
                                 smLocalIsometricEmbedding.cpp
                                 smNonlinearIsometricEmbedding.cpp
                                 smLLE.cpp smIsomap.cpp )
ENDIF(SOLMEC_USE_ARPACK  AND SOLMEC_USE_MKL)

SET(DIMREDUCTION_SOURCES ${DIMREDUCTION_SOURCES_H} ${DIMREDUCTION_SOURCES_CPP})

#------------------------------------------------------------
# LIBRARY smDimensionReduction
#------------------------------------------------------------
ADD_LIBRARY(smDimensionReduction ${DIMREDUCTION_SOURCES})

TARGET_LINK_LIBRARIES(  smDimensionReduction
                        smCommon
                        smIO
                        smUtilities
                        smNumerics
                        smMolecular
                        smLinearSolvers
                        smIterativeSolvers )

IF(SOLMEC_USE_ARPACK)
    TARGET_LINK_LIBRARIES( smDimensionReduction smArpack )
ENDIF(SOLMEC_USE_ARPACK)

IF(UNIX)
    TARGET_LINK_LIBRARIES (smDimensionReduction -lm )
ENDIF(UNIX)

IF(SOLMEC_LIBRARY_PROPERTIES)
  SET_TARGET_PROPERTIES(smDimensionReduction PROPERTIES ${SOLMEC_LIBRARY_PROPERTIES})
ENDIF(SOLMEC_LIBRARY_PROPERTIES)

IF(NOT SOLMEC_INSTALL_NO_LIBRARIES)
  INSTALL(TARGETS smDimensionReduction
    RUNTIME DESTINATION ${SOLMEC_INSTALL_BIN_DIR_CM24} COMPONENT RuntimeLibraries
    LIBRARY DESTINATION ${SOLMEC_INSTALL_LIB_DIR_CM24} COMPONENT RuntimeLibraries
    ARCHIVE DESTINATION ${SOLMEC_INSTALL_LIB_DIR_CM24} COMPONENT Development)
ENDIF(NOT SOLMEC_INSTALL_NO_LIBRARIES)

IF(NOT SOLMEC_INSTALL_NO_DEVELOPMENT)
  FILE(GLOB __files1 "${CMAKE_CURRENT_SOURCE_DIR}/*.h")
  FILE(GLOB __files2 "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
  INSTALL(FILES ${__files1} ${__files2}
    DESTINATION ${SOLMEC_INSTALL_INCLUDE_DIR_CM24}/Numerics/DimensionReduction
    COMPONENT Development)
ENDIF(NOT SOLMEC_INSTALL_NO_DEVELOPMENT)
