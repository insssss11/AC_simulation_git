# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
# ADD_COMPILE_OPTIONS(-Wall -Werror -Wno-shadow=global)
ADD_COMPILE_OPTIONS(-Wall -Wno-shadow=global)
project(AC_simulation)

#----------------------------------------------------------------------------
# Find Geant4, ROOT package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE} ${ROOT_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include
                    ${Geant4_INCLUDE_DIR}
                    ${ROOT_INCLUDE_DIRS})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
set(DESIGN "DEFAULT" CACHE STRING "To changed the detector design, set this variable as TYPE+number.")
set(REFLECTOR "DEFAULT" CACHE STRING "You can select reflector material among ESR, anluminized mylar, tyvek")
set(OPTMODEL "DEFAULT" CACHE STRING "You can select optical boundary model among LUT(or LUTDAVIS), unified")


#----------------------------------------------------------------------------
# geometry options
if(${DESIGN} MATCHES "DEFAULT")
  add_definitions(-DAC_DETECTOR_DESIGN_DEFAULT)

elseif(${DESIGN} MATCHES "LEPS2")
  add_definitions(-DAC_DETECTOR_DESIGN_LEPS2)

elseif(${DESIGN} MATCHES "TYPE20")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE20)
elseif(${DESIGN} MATCHES "TYPE21")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE21)
elseif(${DESIGN} MATCHES "TYPE22")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE22)
elseif(${DESIGN} MATCHES "TYPE23")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE23)
elseif(${DESIGN} MATCHES "TYPE24")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE24)
elseif(${DESIGN} MATCHES "TYPE25")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE25)
elseif(${DESIGN} MATCHES "TYPE26")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE26)
elseif(${DESIGN} MATCHES "TYPE27")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE27)
elseif(${DESIGN} MATCHES "TYPE28")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE28)
elseif(${DESIGN} MATCHES "TYPE29")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE29)
elseif(${DESIGN} MATCHES "TYPE30")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE30)
elseif(${DESIGN} MATCHES "TYPE31")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE31)
elseif(${DESIGN} MATCHES "TYPE32")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE32)
elseif(${DESIGN} MATCHES "TYPE33")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE33)
elseif(${DESIGN} MATCHES "TYPE34")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE34)
elseif(${DESIGN} MATCHES "TYPE35")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE35)
elseif(${DESIGN} MATCHES "TYPE36")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE36)
elseif(${DESIGN} MATCHES "TYPE37")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE37)
elseif(${DESIGN} MATCHES "TYPE38")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE38)
elseif(${DESIGN} MATCHES "TYPE39")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE39)
elseif(${DESIGN} MATCHES "TYPE40")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE40)

elseif(${DESIGN} MATCHES "TYPE1")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE1)
elseif(${DESIGN} MATCHES "TYPE2")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE2)
elseif(${DESIGN} MATCHES "TYPE3")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE3)
elseif(${DESIGN} MATCHES "TYPE4")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE4)
elseif(${DESIGN} MATCHES "TYPE5")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE5)
elseif(${DESIGN} MATCHES "TYPE6")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE6)
elseif(${DESIGN} MATCHES "TYPE7")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE7)
elseif(${DESIGN} MATCHES "TYPE8")
  add_definitions(-DAC_DETECTOR_DESIGN_TYPE8)

elseif(${DESIGN} MATCHES "TESTBENCH1")
  add_definitions(-DAC_TESTBENCH1)
elseif(${DESIGN} MATCHES "TESTBENCH2")
  add_definitions(-DAC_TESTBENCH2)
elseif(${DESIGN} MATCHES "TESTBENCH3")
  add_definitions(-DAC_TESTBENCH3)
elseif(${DESIGN} MATCHES "TESTBENCH4")
  add_definitions(-DAC_TESTBENCH4)

else()
  message(FATAL_ERROR "
-------------------------------------------------------------------------------------
// The variable DESIGN has invalid type number. Option : TYPE1, TYPE2 ..., TYPE8
// Please add an option when running cmake(Default configuration is -DDESIGN=TYPE1)
// (example) cmake .. -DDESIGN=TYPE1
-------------------------------------------------------------------------------------  
")
endif()

#----------------------------------------------------------------------------
# reflector and optical model options
if(${REFLECTOR} MATCHES "DEFAULT")
  add_definitions(-DAC_REFLECTOR_DAVISLUT_ESR)
elseif(${REFLECTOR} MATCHES "DAVISLUT_ESR")
  add_definitions(-DAC_REFLECTOR_DAVISLUT_ESR)
elseif(${REFLECTOR} MATCHES "LUT_TYVEK")
  add_definitions(-DAC_REFLECTOR_LUT_TYVEK)
 
elseif(${REFLECTOR} MATCHES "UNIFIED_ESR")
  add_definitions(-DAC_REFLECTOR_UNIFIED_ESR)
elseif(${REFLECTOR} MATCHES "UNIFIED_ALMYLAR")
  add_definitions(-DAC_REFLECTOR_UNIFIED_ALMYLAR)

else()
  message(FATAL_ERROR "
-------------------------------------------------------------------------------------
// The variable REFLECTOR has invalid string value. Option : ESR and ALMYLAR, TYVEK
// Please add an option when running cmake(Default configuration is -DREFLECTOR=DAVISLUT_ESR)
// (example) cmake .. -DREFLECTOR=DAVISLUT_ESR
-------------------------------------------------------------------------------------  
")
endif()

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDE
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(simulateAC simulateAC.cc ${sources} ${headers})
target_link_libraries(simulateAC ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SIMULATEAC_SCRIPTS
  simulateAC.in
  simulateAC.out
  init_vis.mac
  
  macro/posrun/posrun_main.mac
  macro/posrun/posrun_single.mac
  macro/leps2run/leps2run_main.mac
  macro/tohokurun/tohokurun_main.mac
  macro/testbench/testbench_main.mac
  macro/test_para/test_para_main.mac
  macro/test_para/test_para_single.mac

  vis.mac
  )

foreach(_script ${SIMULATEAC_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(AC_simulation DEPENDS simulateAC)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS simulateAC DESTINATION bin)
message(STATUS "////////////////////////////////////////////////////////////////////////////////////////////")
message(STATUS "// The AC design is ${DESIGN}")

if(${DESIGN} MATCHES "DEFAULT")
  message(STATUS "// Design DEFAULT : DEFAULT design")

elseif(${DESIGN} MATCHES "LEPS2")
  message(STATUS "// Design LEPS2: Full AC construction for LEPS2/SPring8 experiment (TYPE20, 21, 22 used)")  

elseif(${DESIGN} MATCHES "TYPE20")
  message(STATUS "// Design TYPE20: The 30-segment-design to be located at the TPC edge(ESR, TYPEA, DAVIS model)")
elseif(${DESIGN} MATCHES "TYPE21")
  message(STATUS "// Design TYPE21: The 30-segment-design to be located at the TPC edge(ESR, TYPEB, DAVIS model)")
elseif(${DESIGN} MATCHES "TYPE22")
  message(STATUS "// Design TYPE22: The 30-segment-design to be located at the TPC edge(ESR, TYPEC, DAVIS model)")
elseif(${DESIGN} MATCHES "TYPE23")
  message(STATUS "// Design TYPE23: Same with TYPE20 but reflector is Aluminized mylar(TYPEA)")
elseif(${DESIGN} MATCHES "TYPE24")
  message(STATUS "// Design TYPE24: Same with TYPE20 but reflector on Holder is diffusive material(TYPEA)")
elseif(${DESIGN} MATCHES "TYPE25")
  message(STATUS "// Design TYPE25: From type 20, the angle of plane of light guide is 30 degrees, not 45 degrees(TYPEA)")
elseif(${DESIGN} MATCHES "TYPE26")
  message(STATUS "// Design TYPE26: From type 20, the angle of plane of light guide is 0 degrees, not 45 degrees(TYPEA)")
elseif(${DESIGN} MATCHES "TYPE27")
  message(STATUS "// Design TYPE27: From type31, the width is changed to 100 mm")
elseif(${DESIGN} MATCHES "TYPE28")
  message(STATUS "// Design TYPE28: The reflectors on both the holder and light guide are diffuse materials(TYPEA, tyvek, LUT model)")
elseif(${DESIGN} MATCHES "TYPE29")
  message(STATUS "// Design TYPE29: The reflectors on both the holder and light guide are diffuse materials(TYPEB, tyvek, LUT model)")
elseif(${DESIGN} MATCHES "TYPE30")
  message(STATUS "// Design TYPE30: The reflectors on both the holder and light guide are diffuse materials(TYPEC, tyvek, LUT model)")
elseif(${DESIGN} MATCHES "TYPE31")
  message(STATUS "// Design TYPE31: From type 20, Aerogel is not placed inside(TYPEA, ESR, LUT model)")
elseif(${DESIGN} MATCHES "TYPE32")
  message(STATUS "// Design TYPE32: From type 21, Aerogel is not placed inside(TYPEB, ESR, LUT model)")
elseif(${DESIGN} MATCHES "TYPE33")
  message(STATUS "// Design TYPE33: From type 22, Aerogel is not placed inside(TYPEC, ESR, LUT model)")   
elseif(${DESIGN} MATCHES "TYPE34")
  message(STATUS "// Design TYPE34: The 30-segment-design to be located at the TPC edge(ESR, TYPEA, unified model).")
elseif(${DESIGN} MATCHES "TYPE35")
  message(STATUS "// Design TYPE35: The 30-segment-design to be located at the TPC edge(ESR, TYPEB, unified model).") 
elseif(${DESIGN} MATCHES "TYPE36")
  message(STATUS "// Design TYPE36: The 30-segment-design to be located at the TPC edge(ESR, TYPEC, unified model).") 
elseif(${DESIGN} MATCHES "TYPE37")
  message(STATUS "// Design TYPE37: The reflectors on light guide are diffuse materials(TYPEA, tyvek, unified model).") 
elseif(${DESIGN} MATCHES "TYPE38")
  message(STATUS "// Design TYPE38: The reflectors on light guide are diffuse materials(TYPEB, tyvek, unified model).") 
elseif(${DESIGN} MATCHES "TYPE39")
  message(STATUS "// Design TYPE39: The reflectors on light guide are diffuse materials(TYPEC, tyvek, unified model).") 
elseif(${DESIGN} MATCHES "TYPE40")
  message(STATUS "// Design TYPE40: Special type-M with MPPC.")

elseif(${DESIGN} MATCHES "TYPE1")
  message(STATUS "// Design TYPE1 : the original AC detector design")
elseif(${DESIGN} MATCHES "TYPE2")
  message(STATUS "// Design TYPE2 : with the light guide rotated by 180 degrees")
elseif(${DESIGN} MATCHES "TYPE3")
  message(STATUS "// Design TYPE3 : with aerogel alignments shifted")
elseif(${DESIGN} MATCHES "TYPE4")
  message(STATUS "// Design TYPE4 : with aerogel alignments shifted and the light guide rotated by 180 degrees")
elseif(${DESIGN} MATCHES "TYPE5")
  message(STATUS "// Design TYPE5 : with the narrower width")
elseif(${DESIGN} MATCHES "TYPE6")
  message(STATUS "// Design TYPE6 : with the narrower width and the light guide rotated by 180 degrees")
elseif(${DESIGN} MATCHES "TYPE7")
  message(STATUS "// Design TYPE7 : with aerogel alignments shifted and the narrower width")
elseif(${DESIGN} MATCHES "TYPE8")
  message(STATUS "// Design TYPE8 : with aerogel alignments shifted and the narrower width and the light guide rotated by 180 degrees")
elseif(${DESIGN} MATCHES "TESTBENCH1")
  message(STATUS "// Design testbench setup 1 : LED - aerogel - aerogel PMT")
elseif(${DESIGN} MATCHES "TESTBENCH2")
  message(STATUS "// Design testbench setup 2 : LED - aerogel - empty PMT")
elseif(${DESIGN} MATCHES "TESTBENCH3")
  message(STATUS "// Design testbench setup 3 : LED - empty - aerogel PMT")
elseif(${DESIGN} MATCHES "TESTBENCH4")
  message(STATUS "// Design testbench setup 4 : LED - empty - empty PMT")

endif()

message(STATUS "// Reflector tpye and optical model : ${REFLECTOR} ")
message(STATUS "////////////////////////////////////////////////////////////////////////////////////////////")
