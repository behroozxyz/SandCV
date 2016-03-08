#------------------------------------------------------------
# Corral
#------------------------------------------------------------
add_library(Corral Corral.cpp)
target_link_libraries(Corral IO ANNSearch ${MKL_LIBRARIES})

##------------------------------------------------------------
## Multi-Patch SandCV
##------------------------------------------------------------
#add_library(SandCVMultiPatch SandCVMultiPatch.cpp)
#target_link_libraries(SandCVMultiPatch Alignment ClosestPoint)

#------------------------------------------------------------
# Adaptive Biasing Force
#------------------------------------------------------------
add_library(ABF ABF.cpp)
target_link_libraries(ABF Corral SandCV)

#------------------------------------------------------------
# Steering Molecular Dynamics (pulling)
#------------------------------------------------------------
add_library(SMD SMD.cpp)
target_link_libraries(SMD SandCV)
