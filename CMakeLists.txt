cmake_minimum_required( VERSION 3.15 )

project( "Global search algorithm" )

message( STATUS "Global search algorithm" )

if(WIN32)
    option(gtest_force_shared_crt "" TRUE)
endif(WIN32)

include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Configure.cmake" )

############################### GTest ###############################
message( STATUS "############################### GTest ###############################" )
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gtest")
include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest/googletest/include")
if( WIN32 )
    include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest-mpi/win")
else( WIN32 )
    include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest-mpi/linux")
endif( WIN32 )

############################### OpenMP ##############################
message( STATUS "############################### OpenMP ##############################" )
option(USE_OMP ON)
if( USE_OMP OR USE_SEQ )
    find_package( OpenMP )
    if( OpenMP_FOUND )
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else( OpenMP_FOUND )
        set( USE_OMP OFF )
    endif( OpenMP_FOUND )
endif( USE_OMP OR USE_SEQ )

############################ std::thread ############################
message( STATUS "############################ std::thread ############################" )
option(USE_STD ON)
if( USE_STD )
    find_package( Threads )
endif( USE_STD )

############################## CPPLINT ##############################
message( STATUS "############################## CPPLINT ##############################" )
find_package(Python3 REQUIRED)
enable_testing()
if( NOT WIN32 )
    include(cmake/Cpplint.cmake)
    cpplint_add_subdirectory(modules)
endif( NOT WIN32 )

############################## Modules ##############################
message( STATUS "############################## Modules ##############################" )
include_directories(3rdparty/unapproved)
add_subdirectory(modules)