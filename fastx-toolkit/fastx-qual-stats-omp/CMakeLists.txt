set(PRJ_NAME "fastx-qual-stats-omp")
set(LIB_NAME "libfastx")

file(GLOB_RECURSE SOURCES_CPP "${CMAKE_CURRENT_LIST_DIR}/*.cpp")
file(GLOB_RECURSE SOURCES_HPP "${CMAKE_CURRENT_LIST_DIR}/*.hpp")
file(GLOB_RECURSE SOURCES_HXX "${CMAKE_CURRENT_LIST_DIR}/*.hxx")

add_executable(${PRJ_NAME} ${SOURCES_CPP})

find_package(OpenMP REQUIRED)

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PRJ_NAME} PRIVATE OpenMP::OpenMP_CXX)
    target_compile_options(${PRJ_NAME} PRIVATE "/openmp:llvm")
endif()

target_link_libraries(${PRJ_NAME} PRIVATE ${LIB_NAME})
target_include_directories(${PRJ_NAME} PRIVATE "../${LIB_NAME}")