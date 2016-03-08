#------------------------------------------------------------
# SandCV with procrustes superimposition
#------------------------------------------------------------
add_library(SandCV SandCV.cpp)
target_link_libraries(SandCV Alignment ClosestPoint)

#------------------------------------------------------------
# SandCV with contact map
#------------------------------------------------------------
add_library(SandCVCMap SandCVCMap.cpp)
target_link_libraries(SandCVCMap Alignment ClosestPoint)

##------------------------------------------------------------
## Dihedral Angles
##------------------------------------------------------------
#add_library(Dihedral DihedralAngles.cpp)
#target_link_libraries(Dihedral)

#------------------------------------------------------------
# Puckering
#------------------------------------------------------------
add_library(Puckering Puckering.cpp)
target_link_libraries(Puckering ${BOOST_LIBRARIES})
