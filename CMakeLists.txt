cmake_minimum_required(VERSION 3.6.3)

project(mc-pricer CXX)

#############################################################
#                      COMPILATION                          #
#############################################################

# Prefix
if (NOT PREFIX)
        set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
        set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
        message(STATUS "Setting build type to 'Debug' as none was specified.")
        set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()

# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
        message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
        message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()


#############################################################
#                       LIBRARIES                           #
#############################################################

# Detect Boost
set(Boost_NO_BOOST_CMAKE ON)
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
message(STATUS "Boost include directory ${Boost_INCLUDE_DIR}")

# Pnl Directory
## À l'Ensimag
set(Pnl_DIR "/matieres/5MMPMP6/pnl/lib/pnl")
## En local
#set(Pnl_DIR "../pnl/build/lib/pnl")

# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
        add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()


# Detect MPI
set(MPI_LIBS)
find_package(MPI REQUIRED)
if (MPI_CXX_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set (MPI_LIBS ${MPI_LIBS} ${MPI_CXX_LIBRARIES})
else (MPI_CXX_FOUND)
    message(STATUS "MPI not found.")
endif (MPI_CXX_FOUND)

ADD_DEFINITIONS(-DOMPI_SKIP_MPICXX)

#############################################################
#                       SUBDIRECTORIES                      #
#############################################################

add_subdirectory(src)

