set(PRJ_NAME "fastq-to-fasta")
set(LIB_NAME "libfastx")

file(GLOB_RECURSE SOURCES_CPP "${CMAKE_CURRENT_LIST_DIR}/*.cpp")
file(GLOB_RECURSE SOURCES_HPP "${CMAKE_CURRENT_LIST_DIR}/.*hpp")
File(GLOB_RECURSE SOURCES_HXX "${CMAKE_CURRENT_LIST_DIR}/.*hxx")

add_executable(${PRJ_NAME} ${SOURCES_CPP})
target_link_libraries(${PRJ_NAME} PRIVATE ${LIB_NAME})
target_include_directories(${PRJ_NAME} PRIVATE "../${LIB_NAME}")