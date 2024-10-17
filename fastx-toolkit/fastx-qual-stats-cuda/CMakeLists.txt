set(PRJ_NAME "fastx-qual-stats-cuda")
set(LIB_FASTX "libfastx")

set(CXX_STANDARD 20)
set(CUDA_STANDARD 20)
set(CUDA_ARCHITECTURES "89")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
enable_language(CUDA)

file(GLOB_RECURSE CUDA_SOURCES "${CMAKE_CURRENT_LIST_DIR}/*.cu")
file(GLOB_RECURSE SOURCES_CPP "${CMAKE_CURRENT_LIST_DIR}/*.cpp")
file(GLOB_RECURSE SOURCES_HPP "${CMAKE_CURRENT_LIST_DIR}/.*hpp")
File(GLOB_RECURSE SOURCES_HXX "${CMAKE_CURRENT_LIST_DIR}/.*hxx")

add_executable(${PRJ_NAME} ${CUDA_SOURCES})

target_include_directories(${PRJ_NAME} PRIVATE "../${LIB_FASTX}")
target_link_libraries(${PRJ_NAME} PRIVATE ${LIB_FASTX})

target_include_directories(${PRJ_NAME} PRIVATE "../fmt/include")
target_link_libraries(${PRJ_NAME} PRIVATE fmt::fmt)