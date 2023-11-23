cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

list(APPEND CMAKE_PREFIX_PATH "/usr/local/Geant4/geant4-v11.1.3-install/lib/cmake/Geant4")

project(Simulation)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hh)

file(GLOB MACRO_FILES
  "*.mac"
)
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

add_executable(scat scat.cc ${sources} ${headers})

target_link_libraries(scat ${Geant4_LIBRARIES})

add_custom_target(Simulation DEPENDS scat)
