# Project name
set(PROJECT_NAME ${CMAKE_PROJECT_NAME})

# source list
aux_source_directory(. SRC_LIST)
aux_source_directory(./Event-driven_Molecular_Dynamics SRC_LIST)
aux_source_directory(./Forward_Flux_Sampling SRC_LIST)

# library output path
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

# create dynamic library
add_library(EDMD SHARED ${SRC_LIST})
