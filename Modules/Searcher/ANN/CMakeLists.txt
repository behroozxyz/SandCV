#------------------------------------------------------------------------------
# Aproximate Nearest Neighbor Searching (ANN)
#------------------------------------------------------------------------------
# Set the source files of ANN Library
set(ANN_SRCS
    ANN.cpp
    brute.cpp
    kd_tree.cpp
    kd_util.cpp
    kd_split.cpp
    kd_dump.cpp
    kd_search.cpp
    kd_pr_search.cpp
    kd_fix_rad_search.cpp
    bd_tree.cpp
    bd_search.cpp
    bd_pr_search.cpp
    bd_fix_rad_search.cpp
    perf.cpp
    )
# Set the source files
set(ANNSEARCH_SRCS ANNSearch.cpp ${ANN_SRCS})

# Build the library
add_library(ANNSearch ${ANNSEARCH_SRCS})
target_link_libraries(ANNSearch)
