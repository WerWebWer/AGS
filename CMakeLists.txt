cmake_minimum_required (VERSION 2.6) 

get_filename_component(ProjectId ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project (${ProjectId}) 
 
set (HEADERS func.h)
set (SOURCES func.cpp main.cpp) 
 
add_executable (${ProjectId} ${HEADERS} ${SOURCES})