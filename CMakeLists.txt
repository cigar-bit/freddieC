cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED true)
find_package( Threads )
file(GLOB SOURCES CONFIGURE_DEPENDS "*.h" "*.cpp")



add_executable(freddie_segment ${SOURCES})

target_link_libraries(freddie_segment ${CMAKE_THREAD_LIBS_INIT})