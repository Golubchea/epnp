 cmake_minimum_required(VERSION 2.8)

get_filename_component(ProjectId ${CMAKE_CURRENT_LIST_DIR} NAME)
string(REPLACE " " "_" ProjectId ${ProjectId})
project(${ProjectId}  DESCRIPTION "perspective and point problem solved by Lepetit")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)



set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fpermissive -O3 -Wall -Wextra ")
include_directories(  include)
file (GLOB_RECURSE SRC_LIST src/*.cpp src/*.c )

add_library(${PROJECT_NAME} SHARED ${SRC_LIST})
target_link_libraries(${PROJECT_NAME} m)
 
#set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER  ${PROJECT_SOURCE_DIR}/include/.)

