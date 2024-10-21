

# FASTX-Toolkit
This project is based on [agordon/fastx_toolkit](https://github.com/agordon/fastx_toolkit)  
Goal of project is improve processing performance.  
It provides improved performance through 3 ways:

1. `Buffering`
2. `OpenMP`
3. `CUDA`
4. `Cluster`

All tools use buffering by default.  
`OpenMP`, `CUDA`, `Cluster` versions also use buffering by default.
`Cluster` versions based on `OpenMP`.

# Features
- [x] [FASTQ-TO-FASTA](fastx-toolkit/fastq-to-fasta)
- [ ] FASTX Statistics
	- [x] [Buffering](fastx-toolkit/fastx-qual-stats)
	- [x] [OpenMP](fastx-toolkit/fastx-qual-stats-omp)
	- [ ] [CUDA](fastx-toolkit/fastx-qual-stats-cuda)
	- [ ] [Cluster](fastx-toolkit/fastx-qual-stats-cluster)
- [x] [FASTX Sample Generator](fastx-toolkit/fastx-samp-gen)

# Notice
The "Buffering" does not used in C runtime file i/o directly.
Because C runtime file i/o has an internal i/o buffer, the size of which is determined as shown in the following pseudo code:

```Pseudo
/* from glibc */

DiskBlockSize = GetDiskBlockSize()

if DiskBlockSize > 0 && INTERNAL_NO_ERROR
	InternalBufferSize = DiskBlockSize
else
	InternalBufferSize = 1024
```

When you use C standard file i/o functions, like fgets, fread etc... 
A memory copy occurs inside the function, as shown in the following pseudocode:

```
/* from glibc */

any fread(user_defined_buffer, user_defined_buffer_size <- i/obufs)
{
	num_read_bytes = fread_internal(internal_buffer, internal_buffer_size);

	if(num_read_bytes > 0)
		memcpy(user_defined_buffer, internal_buffer, num_read_bytes);
}
```

This means file i/o performance is depend on your disk block size.
So `i/obufs` is not directly applied to C runtime file i/o.

# Usage
You can see help message when you execute program with "-h" flag.  

`FASTQ-TO-FASTA: Convert FASTQ file to FASTA`
|  Option  | Description | Default | Range | 
|:--------:|:-----------:|:-------:|:-----:|
| -h       | print help  |         |       |
| -i       | set input file name | STDIN ||
| -o       | set output file name | STDOUT ||
| -n       | keep sequence with unknown (N) nucleotides.<br/>Default is to discard such sequences. | false ||
| -r       | rename sequence id to number | false ||
| -\-ibufs | set input buffer size | 32768 | IBUFS >= MXSL |
| -\-obufs | set output buffer size | 32768 | > 0 |
| -\-mxsl  | set maximum sequence length | 25000 | > 0 |

`FASTX Sample Generator: Generate FASTX sample`
|  Option  | Description | Default | Range | 
|:--------:|:-----------:|:-------:|:-----:|
| -h        | print help  |         |       |
| -s, -\-sf | set sample format. fasta or fastq |||
| -\-nr     | set num of records || > 0 |
| -\-crlf   | set line break as CRLF | LF ||
| -o        | set output file name | STDOUT ||
| -\-bq     | set base quality offset | 33 ||
| -\-mnq    | set min quality | -15 | BQ + \|MNQ\| >= 0 |
| -\-mxq    | set max quality | 93 | BQ + \|MXQ\| >= MNQ |
| -\-mns    | set min seq length | 1 | > 0 |
| -\-mxs    | set max seq length | 50 | >= MNS |
| -\-obufs  | set output buffer size | 32768 | > 0 |

`FASTX Statistics(Buffering)`
|  Option  | Description | Default | Range | 
|:--------:|:-----------:|:-------:|:-----:|
| -h       | print help  |         |       |
| -i       | set input file name | STDIN ||
| -o       | set output file name | STDOUT ||
| -\-bq    | set base quality offset | 33 | 0 - 255 |
| -\-mnq   | set min quality | -15 | BQ + MNQ >= 0 |
| -\-mxq   | set max quality | 93  | BQ + MXQ <= 255 |
| -\-ibufs | set input buffer size | 32768 | IBUFS >= MXSL |
| -\-mxsl  | set max sequence length | 25000 | > 0 |

`FASTX Statistics(OpenMP)`
|  Option  | Description | Default | Range | 
|:--------:|:-----------:|:-------:|:-----:|
| -h       | print help  |         |       |
| -i       | set input file name | STDIN ||
| -o       | set output file name | STDOUT ||
| -\-bq    | set base quality offset | 33 | 0 - 255 |
| -\-mnq   | set min quality | -15 | BQ + MNQ >= 0 |
| -\-mxq   | set max quality | 93  | BQ + MXQ <= 255 |
| -\-ibufs | set input buffer size | 32768 | IBUFS >= MXSL |
| -\-mxsl  | set max sequence length | 25000 | > 0 |
| -\-rps   | record pool size | 500 ||

# Benchmarks
Device: GA403UI-QS091  
Method: TAKE MINS & ROUND

![Performance Comparison](fastx-toolkit/tests/results/fastx-statistics.png)

| Record Size | Old (ms) | Buffering (ms) | OpenMP (ms) | Buffering Speedup (%) | OpenMP Speedup (%) | OpenMP vs Buffering Increase (%) |
|:-----------:|:--------:|:--------------:|:-----------:|:---------------------:|:-------------------:|:--------------------------------:|
| 1M          | 5259     | 865            | 533         | 83.6                  | 89.9                | 38.4                             |
| 2.5M        | 12974    | 2123           | 1247        | 83.7                  | 90.4                | 41.0                             |
| 5M          | 27106    | 4528           | 2370        | 83.3                  | 91.2                | 47.7                             |
| 10M         | 52969    | 8256           | 4601        | 84.4                  | 91.3                | 44.2                             |
| 30M         | 157420   | 24489          | 13411       | 84.4                  | 91.5                | 45.3                             |
| 50M         | 262089   | 40974          | 21974       | 84.4                  | 91.6                | 46.4                             |
| 100M        | 528681   | 85367          | 52954       | 83.8                  | 89.9                | 38.0                             |