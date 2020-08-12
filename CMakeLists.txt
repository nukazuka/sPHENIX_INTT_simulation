
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(B4c)

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

#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#find_package(ROOT COMPONENTS RooStats)
#include(${ROOT_USE_FILE})

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
#ROOT_GENERATE_DICTIONARY(G__Event Event.h LINKDEF EventLinkDef.h)

#---Create a shared library with geneated dictionary
#add_library(Event SHARED Event.cxx G__Event.cxx)

#---Create  a main program using the library
#add_executable(Main MainEvent.cxx)
#target_link_libraries(Main Event)




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
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable( cal exampleB4c.cc ${sources} ${headers})
#target_link_libraries(cal ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
target_link_libraries(cal ${Geant4_LIBRARIES} )

#----------------------------------------------------------------------------
# Compiler options
# set(CMAKE_CXX_FLAGS "-O1")
 
# When an optimization option is on, the following error is returned:
# [ 11%] Building CXX object CMakeFiles/cal.dir/exampleB4c.cc.o
# In file included from /home/gnukazuk/soft/Geant4/pro/include/Geant4/G4GenericMessenger.hh:35:0,
#                  from /home/gnukazuk/cernbox/Analysis_new/misc/201904_FNAL_test_beam/local/MC/FNAL_calorimeter/include/B4cDetectorConstruction.hh:50,
#                  from /home/gnukazuk/cernbox/Analysis_new/misc/201904_FNAL_test_beam/local/MC/FNAL_calorimeter/exampleB4c.cc:30:
# /home/gnukazuk/soft/Geant4/pro/include/Geant4/G4AnyMethod.hh:74:29: error: ‘remove_const’ is not a class template
#  template<typename T> struct remove_const<const T> {typedef T type;};
#                              ^~~~~~~~~~~~
# CMakeFiles/cal.dir/build.make:62: recipe for target 'CMakeFiles/cal.dir/exampleB4c.cc.o' failed


#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4c. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(EXAMPLEB4C_SCRIPTS
  exampleB4c.out
  exampleB4.in
  gui.mac
  init_vis.mac
#  plotHisto.C
  plotNtuple.C
  run1.mac
#  run2.mac
  run_electron.mac
  vis.mac
  )

foreach(_script ${EXAMPLEB4C_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS cal DESTINATION bin)
