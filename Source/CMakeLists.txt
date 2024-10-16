add_library(${NSEOF_PROJECT_NAME} OBJECT)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*")
list(FILTER SOURCES EXCLUDE REGEX ".*Main\\.cpp$")

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})

target_sources(${NSEOF_PROJECT_NAME} PRIVATE ${SOURCES})

target_link_libraries(${NSEOF_PROJECT_NAME} PUBLIC NS-EOF-Interface)
target_include_directories(${NSEOF_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(${NSEOF_PROJECT_NAME}-Runner Main.cpp)
target_link_libraries(${NSEOF_PROJECT_NAME}-Runner PRIVATE ${NSEOF_PROJECT_NAME})

if(MSVC)
  configure_file(${CMAKE_SOURCE_DIR}/petsc_commandline_arg ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
else()
  configure_file(${CMAKE_SOURCE_DIR}/petsc_commandline_arg ${CMAKE_BINARY_DIR} COPYONLY)
endif()
