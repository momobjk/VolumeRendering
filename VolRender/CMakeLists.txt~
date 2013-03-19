project(VolRender)
set(QT_QMAKE_EXECUTABLE /usr/bin/qmake-qt4)
cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)

find_package(Qt4 REQUIRED)
find_package(GLEW REQUIRED)
find_package(GLUT REQUIRED)


#set( CMAKE_VERBOSE_MAKEFILE on )
#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_BUILD_TYPE Debug)
#modules complémentaires (glew et autres)


#LIBRARIES
find_library(X11    X11)
if (X11)
    message("-- Found libX11    (${X11})")
endif (X11)

include_directories(${QT_INCLUDES} ${CMAKE_CURRENT_BINARY_DIR} ${GLUT_INCLUDE_DIR} ${GLEW_INCLUDE_PATH})
include_directories(headers headers/QT_MOC)
include(${QT_USE_FILE})


# OPENMP
#if(OPENMP_FOUND)
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

#c++11 and optimization
if(CMAKE_COMPILER_IS_GNUCXX)
   execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
   if (GCC_VERSION VERSION_GREATER 4.7 OR GCC_VERSION VERSION_EQUAL 4.7)
        message(STATUS "C++11 activated.")
        add_definitions("-std=c++11")
   elseif(GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
        message(WARNING "C++0x activated. If you get any errors update to a compiler which fully supports C++11")
        add_definitions("-std=c++0x")
   else ()
        message(FATAL_ERROR "C++11 needed. Therefore a gcc compiler with a version higher than 4.3 is needed.")   
   endif()
else(CMAKE_COMPILER_IS_GNUCXX)
   add_definitions("-std=c++0x") 
endif(CMAKE_COMPILER_IS_GNUCXX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -march=native")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O2")
set(CMAKE_CXX_FLAGS_RELEASE  "${CMAKE_CXX_FLAGS_RELEASE}")


message("BUILD TYPE : ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  message("COMPILER FLAGS : ${COMPILE_DEFINITIONS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else()
  message("COMPILER FLAGS : ${COMPILE_DEFINITIONS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif()

#The source files
file(GLOB SRC src/*.cpp)
#set(SRCS inpainting.cpp main.cpp view.cpp)
#moc headers
file(GLOB MOC_HDRS headers/QT_MOC/*.h)

#uis
#set(UIS view.ui)
file(GLOB UIS uis/*.ui)


qt4_wrap_cpp(MOC_SRC ${MOC_HDRS}) 

qt4_wrap_ui(UI_HDRS ${UIS})

add_executable(VolRender ${SRC} 
			${MOC_SRC}
			${UI_HDRS})
					

target_link_libraries(VolRender ${QT_LIBRARIES} ${X11} ${GLUT_LIBRARIES} ${GLEW_LIBRARY})
