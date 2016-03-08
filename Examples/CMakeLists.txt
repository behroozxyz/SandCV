#------------------------------------------------------------
# Procrustes Superimposition
#------------------------------------------------------------
add_executable(ProcrustesSuperimpositionCheck ProcrustesSuperimpositionCheck.cpp)
target_link_libraries(ProcrustesSuperimpositionCheck Alignment IO)

#------------------------------------------------------------
# Contact Map
#------------------------------------------------------------
add_executable(ContactMapCheck ContactMapCheck.cpp)
target_link_libraries(ContactMapCheck Alignment IO)

#------------------------------------------------------------
# ANN Searcher
#------------------------------------------------------------
add_executable(ANNSearchCheck ANNSearchCheck.cpp)
target_link_libraries(ANNSearchCheck ANNSearch ${MATLAB_LIBRARIES})

#------------------------------------------------------------
# Local Maximum Entropy Approximant
#------------------------------------------------------------
add_executable(LocalMaxEntCheck LocalMaxEntCheck.cpp)
target_link_libraries(LocalMaxEntCheck LocalMaxEnt)

#------------------------------------------------------------
# Closest Point Projection
#------------------------------------------------------------
add_executable(ClosestPointProjectionCheck ClosestPointProjectionCheck.cpp)
target_link_libraries(ClosestPointProjectionCheck ClosestPoint ${MATLAB_LIBRARIES})

#------------------------------------------------------------
# SandCV
#------------------------------------------------------------
add_executable(SandCVCheck SandCVCheck.cpp)
target_link_libraries(SandCVCheck SandCV IO)

#------------------------------------------------------------
# SandCV with smooth contact map
#------------------------------------------------------------
add_executable(SandCVCMapCheck SandCVCMapCheck.cpp)
target_link_libraries(SandCVCMapCheck SandCVCMap IO)

#------------------------------------------------------------
# Corral
#------------------------------------------------------------
add_executable(CorralCheck CorralCheck.cpp)
target_link_libraries(CorralCheck Corral)

#------------------------------------------------------------
# Puckering Coordinate
#------------------------------------------------------------
add_executable(PuckeringCheck PuckeringCheck.cpp)
target_link_libraries(PuckeringCheck Puckering IO)

#------------------------------------------------------------
# ABF
#------------------------------------------------------------
if(USE_SWIG)
    add_executable(ABFCheck ABFCheck.cpp)
    target_link_libraries(ABFCheck ABF SandABF)
endif(USE_SWIG)
