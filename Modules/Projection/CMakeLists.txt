#------------------------------------------------------------
# Closest Point Projection
#------------------------------------------------------------
add_library(ClosestPoint ClosestPointProjection.cpp)
target_link_libraries(ClosestPoint LocalMaxEnt ANNSearch IO ${MKL_LIBRARIES})

