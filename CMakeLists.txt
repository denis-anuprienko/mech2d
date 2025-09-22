cmake_minimum_required(VERSION 2.8)

project(mech2d)

option(USE_MPI "Compile with MPI support" OFF)

add_executable(mech2d mech2d.cpp)

#set(INMOST_DIR C:/Users/denis/Documents/repos/INMOST/)
find_package(inmost REQUIRED)
if(NOT inmost_FOUND)
    message("INMOST not found!")
endif()

include_directories(${INMOST_INCLUDE_DIRS})
message(${INMOST_INCLUDE_DIRS})
message(${INMOST_LIBRARIES})
target_link_libraries(mech2d ${INMOST_LIBRARIES})
