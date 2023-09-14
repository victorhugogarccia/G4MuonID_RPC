cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(MuonID)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Add include directories
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR})

add_executable(MuonID MuonID_sim.cc ${sources} ${headers})
target_link_libraries(MuonID ${Geant4_LIBRARIES})
