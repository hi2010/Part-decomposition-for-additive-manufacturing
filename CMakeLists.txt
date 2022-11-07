cmake_minimum_required (VERSION 3.2)

project(testPro)
# this build only works with linux based systems because the linked libs are .so and build here and the install is not defined as to copy the libraries where they need to go (depening on the os)
# as a reminder because there is some os checking in this file ...

# on install the files go to top level /bin folder
# if the install fails with cmake/install.util ... failed to open for writing (permission denied)
# it does not matter
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/bin)
# i'm not sure if there rpath things are needed or even work ...
# well they don't work completely ... not sure if at all but it runs crntly so ...
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#region import all needed for lib3mf(core) (CopySharedLibrary call on the target(s) is still needed afterwards)
#set(CMAKE_CXX_STANDARD 11)
set(EXECUTABLE_NAME decomposer)

set(VHACD_DIR ${CMAKE_SOURCE_DIR}/external/VHACD_Performance/VHACD_Lib)
#set(VHACD_DIR ${CMAKE_SOURCE_DIR}/external/VHACD_Performance_1/src/VHACD_Lib)

# determine the lib ending used by the sys (copied from lib3mf) and used to copy lib3mf
# cmake cxx flag get overwritten later in the program
if (WIN32)
	set(LSUFFIX "dll")
	set(LSUFFIXOUT ".dll")
elseif(UNIX AND NOT APPLE)
	set(LSUFFIX "so")
	set(LSUFFIXOUT ".so.2")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2")
else()
	set(LSUFFIX "dylib")
	set(LSUFFIXOUT ".2.dylib")
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2")
endif()

#add_definitions( -DTEXTURESPATH="${CMAKE_CURRENT_SOURCE_DIR}/../Files/Textures/" )

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/lib3mf/Bindings/Cpp)

#if (WIN32)
#	link_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/lib3mf/Lib)
#	link_libraries(lib3mf)
#else()
	# Unix prefixes the name of the library with "lib" anyway
#	link_directories(${CMAKE_CURRENT_SOURCE_DIR}/external/lib3mf/Bin)
#	link_libraries(3mf)
#	link_directories(${CMAKE_CURRENT_BINARY_DIR})
#endif()

function(CopySharedLibrary MyTarget)
	add_custom_command(TARGET ${MyTarget} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/external/lib3mf/Bin/lib3mf.${LSUFFIX}" $<TARGET_FILE_DIR:${MyTarget}>/lib3mf${LSUFFIXOUT})
endfunction(CopySharedLibrary)

if (${MSVC})
	IF(${CMAKE_VERSION} VERSION_LESS 3.6.3)
		MESSAGE ("Note: You need to manually select a StartUp-project in Visual Studio.")
	ELSE()
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Example_Cube)
	ENDIF()
endif()


#include_directories(${VHACD_DIR}/public)
#include_directories(${VHACD_DIR}/inc)
#link_directories(${VHCAD_DIR}/build)
#link_libraries(vhacd)
#link_directories(${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(${VHACD_DIR})

#project(VHACD_TEST)
# include(${CMAKE_COMMON_INC})

if (NOT NO_OPENCL)
	find_package(OpenCL QUIET)
endif()
if (NOT NO_OPENMP)
	find_package(OpenMP QUIET)
endif()
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#endregion

#region import for vtk
# uncomment for newer version glibc2.3 ++
#set(VTK_DIR ${CMAKE_SOURCE_DIR}/external/vtk/build)
#set(VTK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/vtk_old_build/build)
set(VTK_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/vtk/build)
#find_package(VTK REQUIRED)
find_package(VTK COMPONENTS 
  vtkCommonColor
  vtkCommonCore
  vtkCommonTransforms
  vtkIOGeometry
  vtkCommonDataModel
  vtkFiltersFlowPaths  # why not work
  vtkFiltersGeneral
  vtkFiltersSources
  vtkFiltersExtraction
  vtkFiltersModeling
  vtkFiltersGeometry
  vtkIOInfovis
  vtkInteractionStyle
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  QUIET
)
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")
#add_executable(DistanceBetweenPoints MACOSX_BUNDLE DistanceBetweenPoints.cxx )

#endregion

# prbly should have used current source dir everywhere
set(OPENGA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/openGA-master/src)

add_executable(${EXECUTABLE_NAME} ${CMAKE_SOURCE_DIR}/src/main.cpp ${CMAKE_SOURCE_DIR}/src/oclHelper.cpp ${CMAKE_SOURCE_DIR}/src/lVhacdHelper.cpp ${CMAKE_SOURCE_DIR}/src/l3mfHelper.cpp  ${CMAKE_SOURCE_DIR}/src/lVtkHelper.cpp ${CMAKE_SOURCE_DIR}/src/polyQualityMeasures.cpp ${CMAKE_SOURCE_DIR}/src/workflows.cpp ${CMAKE_SOURCE_DIR}/src/modifiedClipClosedSurface.cpp)
# add lib3mf.so.2
CopySharedLibrary(${EXECUTABLE_NAME})
#add_custom_command(TARGET testex POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${VHACD_DIR}/build/libvhacd.a" $<TARGET_FILE_DIR:testex>/libvhacd.a)
# link vhacd
target_include_directories(${EXECUTABLE_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/inc ${VHACD_DIR}/public ${VHACD_DIR}/inc ${OPENGA_DIR})
#target_include_directories(testex PRIVATE ${VHACD_DIR}/inc ${VHACD_DIR}/public)
#link_directories(${VHACD_DIR}/build)
#link_directories(/home/louis/Code/testProj0/external/VHACD_Performance/VHACD_Lib/build)
#target_link_libraries(testex libvhacd)
#find_library(testex libvhacd ${VHACD_DIR}/build)
#target_link_libraries(testex vhacd)
target_link_libraries(${EXECUTABLE_NAME} ${VHACD_DIR}/build/libvhacd.a)
#target_link_libraries(${EXECUTABLE_NAME} 3mf)

list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_SOURCE_DIR}/external/lib3mf/Bin)
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_CURRENT_BINARY_DIR})
find_library(3MFLIB 3mf)
#message("lib3mf: ${3MFLIB}")
target_link_libraries(${EXECUTABLE_NAME} ${3MFLIB})
#target_link_libraries(${EXECUTABLE_NAME} ${CMAKE_CURRENT_BINARY_DIR}/lib3mf.so.2)

#add_library(libvhacd STATIC IMPORTED) # or STATIC instead of SHARED
#set_target_properties(libvhacd PROPERTIES
#  IMPORTED_LOCATION "${VHACD_DIR}/build/libvhacd.a"
#  INTERFACE_INCLUDE_DIRECTORIES "${VHACD_DIR}/public"
#)
#target_link_libraries(testex libvhacd)

#add_executable(testVHACD ${PROJECT_CPP_FILES} ${PROJECT_C_FILES} ${PROJECT_INC_FILES} ${PROJECT_INL_FILES})

#target_link_libraries(testVHACD vhacd)
#target_link_libraries(testex ${VHACD_DIR}/build/libvhacd.a)



#region imports / definitions for vhacd
if (OpenCL_FOUND)
    include_directories("${OpenCL_INCLUDE_DIRS}")
    add_definitions( -DOPENCL_FOUND=1 )
    #target_link_libraries(testVHACD ${OpenCL_LIBRARIES})
    target_link_libraries(${EXECUTABLE_NAME} ${OpenCL_LIBRARIES})
endif()

if (NOT WIN32 AND NOT APPLE)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall")
   #target_link_libraries(testVHACD rt)
   target_link_libraries(${EXECUTABLE_NAME} rt)
endif()
#endregion


# vhacd

#file(GLOB_RECURSE VHACD_SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "/external/VHACD_Performance/VHACD_Lib/src/*.c" "/external/VHACD_Performance/VHACD_Libsrc/*.cpp")
#file(GLOB_RECURSE VHACD_INCLS RELATIVE ${VHACD_DIR} "public/*.h" "public/*.hpp" "inc/*.h" "inc/*.hpp")
#target_include_directories(libvh ${VHACD_DIR})
#find_package()
#MESSAGE (WARNING "srcs ${VHACD_DIR}")
#MESSAGE (WARNING "srcs ${VHACD_SOURCES}")
#add_library(libvhawcd STATIC ${VHACD_SOURCES} ${VHACD_INCLS})
#add_subdirectory(${VHACD_DIR})
#target_link_libraries(testex VHACD_LIB)

target_link_libraries(${EXECUTABLE_NAME} ${VTK_LIBRARIES})
# vtk_module_autoinit is needed
vtk_module_autoinit(
    TARGETS ${EXECUTABLE_NAME}
    MODULES ${VTK_LIBRARIES}
    )

set_target_properties(${EXECUTABLE_NAME} PROPERTIES CXX_STANDARD 17)
# for no cmake build type just use defaults
# change to O0 for debug
if (CMAKE_BUILD_TYPE EQUAL "DEBUG" OR CMAKE_BUILD_TYPE EQUAL "Debug")
    SET(GCC_COVERAGE_COMPILE_FLAGS "-g -pthread -O0")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
endif()
if (CMAKE_BUILD_TYPE EQUAL "RELEASE" OR CMAKE_BUILD_TYPE EQUAL "Release")
    SET(GCC_COVERAGE_COMPILE_FLAGS "-g -pthread -O2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} -O2")
endif()
set(CMAKE_CXX_COMPILER_ID "g++")
#set(CMAKE_CXX_COMPILER_VERSION "9.4")

# so this is a little bit a mess but well idk r8now as the whole cmkae file is a mess
# as a reminder (2 myself):
# anyway ... this copies all files and libs needed to the CMAKE_INSTLL_PREFIX (./bin) directory
# but as i was unable to set the RPATH correct or so use the start script in ./bin folder to run the program if libs can not be found
# which is the case if they're not installed additionally ... so the scripts adds the ./bin/lib folder to the LD_LIRARY_PATH so they can be found
# alternatively one could install them using apt or apk or yum or ...
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath -Wl,./ ../external/vtk/build/lib/${VTK_LIBRARIES}")
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
install(TARGETS ${EXECUTABLE_NAME} DESTINATION bin)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external/vtk/build/lib DESTINATION .)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/build/lib3mf${LSUFFIXOUT} DESTINATION lib)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/build/lib3mf${LSUFFIXOUT} DESTINATION bin)
#install(TARGETS ${VTK_LIBRARIES} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin)

