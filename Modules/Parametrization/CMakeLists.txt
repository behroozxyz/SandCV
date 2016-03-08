#------------------------------------------------------------
# LocalMaxEnt
#------------------------------------------------------------
add_library(LocalMaxEnt LocalMaxEntropy.cpp)
target_link_libraries(LocalMaxEnt ANNSearch ${MKL_LIBRARIES})
