
# FASTX-Toolkit
The project is based on [agordon/fastx_toolkit](https://github.com/agordon/fastx_toolkit)

The goal is to improve processing performance.
It provides improved performance through 3 methods:

1. Buffering
2. OpenMP
3. CUDA 

OpenMP and CUDA versions use buffering by default.

# Features
- [x] FASTX Statistics
	- [x] [Buffering](fastx-toolkit/fastx-qual-stats)
	- [x] [OpenMP](fastx-toolkit/fastx-qual-stats-omp)
	- [ ] [CUDA](fastx-toolkit/fastx-qual-stats-cuda) #1
- [x] [FASTX Sample Generator](fastx-samp-gen)

# Benchmarks
