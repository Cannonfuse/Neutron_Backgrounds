#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
project(SIMPDETPART)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/)

# set(TETGEN_INCLUDE_DIR "/usr/include/")
# set(TETGEN_LIBRARY "/usr/lib/x86_64-linux-gnu/")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(CMAKE_CXX_STANDARD 11)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

# find_package(tetgen REQUIRED)

#----------------------------------------------------------------------------
# Compilier flags
#
include(CMakeDependentOption)

option(ONEINCHCLYC "Use the 1 inch crystal" ON)
if(ONEINCHCLYC)
  add_definitions(-DONEINCHCLYC) 
endif()

option(THREEINCHCLYC "Use the 3 inch crystal" OFF)
if(THREEINCHCLYC)
  add_definitions(-DTHREEINCHCLYC) 
endif()

option(USE_TUNNEL "Use the tunnel structure" OFF)
if(USE_TUNNEL)
  add_definitions(-DUSE_TUNNEL) 
endif()

option(USE_CASING "Use the detector casing" OFF)
if(USE_CASING)
  add_definitions(-DUSE_CASING) 
endif()

option(USE_FTC "Use the fine tunnel collimator" OFF)
if(USE_FTC)
  add_definitions(-DUSE_FTC) 
endif()

option(USE_MTC "Use the medium tunnel collimator" OFF)
if(USE_MTC)
  add_definitions(-DUSE_MTC) 
endif()

option(USE_LTC "Use the large tunnel collimator" OFF)
if(USE_LTC)
  add_definitions(-DUSE_LTC) 
endif()

cmake_dependent_option(USE_SHIELD1 "Use the parafin shield type 1" OFF "NOT USE_SHIELD2" OFF)
if(USE_SHIELD1)
  add_definitions(-DUSE_SHIELD1) 
endif()

cmake_dependent_option(USE_SHIELD2 "Use the parafin shield type 2" OFF "NOT USE_SHIELD1" OFF)
if(USE_SHIELD2)
  add_definitions(-DUSE_SHIELD2) 
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/source/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)

#----------------------------------------------------------------------------
# Locate libraries
# 
#
add_library(OUG4Materials SHARED IMPORTED)
set_target_properties(OUG4Materials PROPERTIES
  IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libOUG4Materials.so"
)
add_library(tetgen SHARED IMPORTED)
set_target_properties(tetgen PROPERTIES
  IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/lib/libtet.so"
)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(neutron_backgrounds neutron_backgrounds.cpp ${sources} ${headers})
target_link_libraries(neutron_backgrounds PUBLIC ${Geant4_LIBRARIES} OUG4Materials tetgen)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB1_SCRIPTS
  exampleSIMPDETPART.in
  exampleSIMPDETPART.out
  init_vis.mac
  run1.mac
  run2.mac
  runE.mac
  vis.mac
  Be9dn.mac
  Uniform.mac
  UniformC7LYC.mac
  Uniform_C6LYC.mac
  Uniform_C7LYC.mac
  Uniform_Be9dn_12m.mac
  Uniform_Be9pn_6m.mac
  Uniform_Ddn_6m.mac
  UniformE_C7LYC_MS_Test.mac
  LowE.mac
  SingleE.mac
  eKin
  run_dn12m_uniform_energies
  run_pn6m_uniform_energies
  )

foreach(_script ${EXAMPLEB1_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

file(COPY CADModels DESTINATION ${PROJECT_BINARY_DIR})
file(COPY lib DESTINATION ${PROJECT_BINARY_DIR})


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(SIMPDETPART DEPENDS neutron_backgrounds)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS neutron_backgrounds DESTINATION bin)
