﻿add_subdirectory(libfastx)
add_subdirectory(fmt)
add_subdirectory(fastq-to-fasta)
add_subdirectory(fastx-qual-stats)
add_subdirectory(fastx-qual-stats-omp)
# add_subdirectory(fastx-qual-stats-cuda)
add_subdirectory(fastx-samp-gen)

if (CMAKE_VERSION VERSION_GREATER 3.20)
	set_property(TARGET libfastx PROPERTY CXX_STANDARD 20)
	set_property(TARGET fmt PROPERTY CXX_STANDARD 20)
	set_property(TARGET fastq-to-fasta PROPERTY CXX_STANDARD 20)
	set_property(TARGET fastx-qual-stats PROPERTY CXX_STANDARD 20)
	set_property(TARGET fastx-qual-stats-omp PROPERTY CXX_STANDARD 20)
	set_property(TARGET fastx-samp-gen PROPERTY CXX_STANDARD 20)
endif()
