cmake_minimum_required (VERSION 2.6)
project (gossip)
include_directories ("${PROJECT_SOURCE_DIR}/src")

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

set(CMAKE_BUILD_TYPE Release)
add_definitions(-DRETIN_DISABLE_JAVA)
add_definitions(-DTHREAD_POOL_USE_POSIX_THREADS)

IF(CMAKE_COMPILER_IS_GNUCC)
	SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}  -g -Wall -std=c99   -pedantic-errors -Wno-long-long -mfpmath=sse -msse")
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-sign-compare -std=c++0x -Wno-long-long -mfpmath=sse -msse")
ENDIF(CMAKE_COMPILER_IS_GNUCC)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib/")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
#SET(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS}   -march=native")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

set(core_dir src/retin/toolbox/core)
file(GLOB core_headers ${core_dir}/core.h ${core_dir}/arrays.h ${core_dir}/auto_array_ptr.h ${core_dir}/string.h ${core_dir}/fvec_reader.h  ${core_dir}/siftgeo_reader.h ${core_dir}/fvec_writer.h ${core_dir}/work_task.h ${core_dir}/thread_pool.h ${core_dir}/work_set.h ${core_dir}/work_queue.h ${core_dir}/work_runnable.h ${core_dir}/SharedMatrix.h ${core_dir}/sharedmem/shared_mem.h ${core_dir}/sharedmem/matrix_writer.h ${core_dir}/sharedmem/matrix_reader.h)
file(GLOB core_sources ${core_dir}/SharedMatrix.cpp ${core_dir}/SharedMatrixDouble.cpp ${core_dir}/sharedmem/matrix_reader.cpp ${core_dir}/sharedmem/matrix_writer.cpp ${core_dir}/sharedmem/shared_mem.cpp)   

set(algebra_dir src/retin/toolbox/algebra)
file(GLOB algebra_headers ${algebra_dir}/core.h ${algebra_dir}/vector_uchar.h ${algebra_dir}/vector_short.h ${algebra_dir}/vector_float.h ${algebra_dir}/vector_double.h ${algebra_dir}/matrix_float.h  ${algebra_dir}/matrix_double.h ${algebra_dir}/vector.h ${algebra_dir}/matrix.h  ${algebra_dir}/pca.h)
file(GLOB algebra_sources ${algebra_dir}/core.c ${algebra_dir}/vector_uchar.c ${algebra_dir}/vector_short.c ${algebra_dir}/vector_float.c ${algebra_dir}/vector_double.c ${algebra_dir}/matrix_float.c ${algebra_dir}/matrix_double.c ${algebra_dir}/pca.cpp )

file(GLOB common_sources src/common/math.cpp src/common/matrix.cpp src/common/multithread.cpp src/common/utils.cpp src/common/matrix_float.cpp src/common/matrix_double.cpp src/common/matrixdouble.cpp)
file(GLOB common_headers src/common/math.h src/common/matrix.h src/common/matrixdouble.h)


add_library(retin_algebra ${algebra_sources})
add_library(retin_core ${core_sources} ${document_sources})

add_executable(gossip_pca src/main.cpp ${common_sources})

set(retin_libs pthread retin_core retin_algebra)

target_link_libraries(gossip_pca ${retin_libs})


find_package(PNG)
if(PNG_FOUND)
	include_directories(${PNG_INCLUDE_DIR})
	target_link_libraries(gossip_pca ${PNG_LIBRARIES} )
	add_definitions(-DRETIN_ENABLE_PNG)
endif(PNG_FOUND)

find_package(JPEG)
if(JPEG_FOUND)
	include_directories(${JPEG_INCLUDE_DIR})
	target_link_libraries(gossip_pca ${JPEG_LIBRARIES} )
	add_definitions(-DRETIN_ENABLE_JPEG)
endif(JPEG_FOUND)

