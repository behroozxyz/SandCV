#------------------------------------------------------------------------------
# Procrustes Superimpostison (ProcSup)
#------------------------------------------------------------------------------
add_library(Alignment ProcrustesSuperimposition.cpp 
                      ContactMap.cpp)
target_link_libraries(Alignment ${MKL_LIBRARIES})
