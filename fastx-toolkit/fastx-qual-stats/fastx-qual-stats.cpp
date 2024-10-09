#include <cinttypes>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <format>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <vector>
#include "args.hxx"
#include "fastx.hpp"

using namespace std;

enum class Nucleotide : uint8_t
{
	ALL, A, C, G, T, N,
	_NUCLEOTIDE_COUNT_,
	UNDEFINED
};

enum class OutputVersion : uint8_t
{
	V1, V2, UNDEFINED,
};

static constexpr uint8_t ALL = static_cast<uint8_t>(Nucleotide::ALL);
static const vector<char> NUC_CHARS = { '\0', 'A', 'C', 'G', 'T', 'N' };
static const vector<string> COMMON_HEADERS = { "count", "min", "max", "sum", "mean", "Q1", "med", "Q3", "IQR", "lW", "rW" };

struct NucleotideStatistics
{
	int min = 100;
	int max = -100;
	uint64_t sum = 0;
	uint64_t count = 0;
	uint64_t base_counts[MAX_QUALITY - MIN_QUALITY] = { 0 };
};

struct ColumnStatistics
{
	NucleotideStatistics nuc_stats[static_cast<uint8_t>(Nucleotide::_NUCLEOTIDE_COUNT_)];
};

/* Argument parser constants */
static const char* PROLOGUE = "";
static const char* EPILOGUE = "";

/* I/O variables */
static char* in_buf;

/* Argument variables */
static OutputVersion out_ver = OutputVersion::V1;

static char base_qual_offset = BASE_QUALITY_OFFSET;
static char min_qual = MIN_QUALITY;
static char max_qual = MAX_QUALITY;

static size_t in_buf_size = 32768;

/* Statistics variables */
static ColumnStatistics* col_stats;
static int nuc_idxs[numeric_limits<uint8_t>::max() + 1] = { 0 };

/* libfastx variables */
static FastxContext_t fastx_ctx;
static FastxRecord_t fastx_record;

void valid_args()
{
	bool result = true;

	result &= min_qual < max_qual;
	result &= fastx_ctx.max_seq_len > 0;
	result &= in_buf_size > 0 && in_buf_size >= fastx_ctx.max_seq_len;

	if (!result)
		throw invalid_argument("Invalid arguments");
}

void parse_args(int argc, char** argv)
{
	args::ArgumentParser parser(PROLOGUE, EPILOGUE);
	args::HelpFlag help(parser, "help", "Display options", { 'h', "help" });

	args::MapFlag<string, OutputVersion> ov_arg(parser, "ov", "output version: v1 or v2", { "ov" }, {
	{ "v1", OutputVersion::V1 },
	{ "v2", OutputVersion::V2 } }, out_ver);

	args::ValueFlag<string> in_arg(parser, "in", "input file name. default is STDIN", { 'i' }, fastx_ctx.in_name);
	args::ValueFlag<string> out_arg(parser, "out", "output file name. default is STDOUT", { 'o' }, fastx_ctx.out_name);

	args::Group qual_group(parser, "Quality");
	args::ValueFlag<char> bq_arg(qual_group, "bq", format("base quality offset. default is {}", static_cast<int>(base_qual_offset)), { "bq" }, base_qual_offset);
	args::ValueFlag<char> mnq_arg(qual_group, "mnq", format("min quality. default is {}", static_cast<int>(min_qual)), { "mnq" }, min_qual);
	args::ValueFlag<char> mxq_arg(qual_group, "mxq", format("max quality. default is {}", static_cast<int>(max_qual)), { "mxq" }, max_qual);

	args::Group io_tuning_group(parser, "I/O Tuning");
	args::ValueFlag<size_t> ibufs_arg(io_tuning_group, "ibufs", format("input buffer size. default is {}", in_buf_size), { "ibufs" }, in_buf_size);
	args::ValueFlag<size_t> mxsl_arg(io_tuning_group, "mxsl", format("max sequence length. default is {}", fastx_ctx.max_seq_len), { "mxsl" }, fastx_ctx.max_seq_len);

	try
	{
		parser.ParseCLI(argc, argv);

		out_ver = args::get(ov_arg);

		args::get(in_arg).copy(fastx_ctx.in_name, MAX_PATH, 0);
		args::get(out_arg).copy(fastx_ctx.out_name, MAX_PATH, 0);

		base_qual_offset = args::get(bq_arg);
		min_qual = args::get(mnq_arg);
		max_qual = args::get(mxq_arg);

		in_buf_size = args::get(ibufs_arg);
		fastx_ctx.max_seq_len = args::get(mxsl_arg);

		valid_args();
	}
	catch (const exception& e)
	{
		cout << parser << endl;
		cerr << e.what() << endl;
		exit(errno);
	}
}

void init_vals()
{
	for (int i = 0; i < NUC_CHARS.size(); ++i)
	{
		nuc_idxs[NUC_CHARS[i]] = i;
		nuc_idxs[tolower(NUC_CHARS[i])] = i;
	}
}

void alloc_bufs()
{
	in_buf = new char[in_buf_size];

	memset(&fastx_record, 0, sizeof(FastxRecord_t));
	fastx_record.seq_id = new char[fastx_ctx.max_seq_len];
	fastx_record.seq = new char[fastx_ctx.max_seq_len];
	fastx_record.qual = new char[fastx_ctx.max_seq_len];

	col_stats = new ColumnStatistics[fastx_ctx.max_seq_len];
}

void free_bufs()
{
	if (in_buf) delete[] in_buf;

	if (fastx_record.seq_id) delete[] fastx_record.seq_id;
	if (fastx_record.seq) delete[] fastx_record.seq;
	if (fastx_record.qual) delete[] fastx_record.qual;

	if (col_stats) delete[] col_stats;
}

void open_files()
{
	open_file(fastx_ctx.in_name, "rb", &fastx_ctx.in_stream);
	fastx_ctx.format = get_file_format(fastx_ctx.in_stream);

	open_file(fastx_ctx.out_name, "wb", &fastx_ctx.out_stream);
}

void close_files()
{
	close_file(fastx_ctx.in_stream);
	close_file(fastx_ctx.out_stream);
}

void update_nuc_statistics(size_t col_idx, uint8_t nuc_idx, int qual, size_t read_count)
{
	col_stats[col_idx].nuc_stats[nuc_idx].count += read_count;

	if (fastx_ctx.format == FileFormat::FILE_FORMAT_FASTQ)
	{
		col_stats[col_idx].nuc_stats[nuc_idx].min = min(col_stats[col_idx].nuc_stats[nuc_idx].min, qual);
		col_stats[col_idx].nuc_stats[nuc_idx].max = max(col_stats[col_idx].nuc_stats[nuc_idx].max, qual);
		col_stats[col_idx].nuc_stats[nuc_idx].sum += qual;
		col_stats[col_idx].nuc_stats[nuc_idx].base_counts[qual - MIN_QUALITY] += read_count;
	}
}

void process_record(FastxRecord_t* record, size_t record_idx)
{
	for (size_t i = 0; i < record->seq_len; ++i)
	{
		char nuc = fastx_record.seq[i];
		int qual = fastx_record.qual[i] - base_qual_offset;

		update_nuc_statistics(i, ALL, qual, record->read_count);
		update_nuc_statistics(i, nuc_idxs[nuc], qual, record->read_count);
	}
}

void read_records()
{
	size_t record_idx = 0; // FIXED: Receive only one record at each time.
	size_t num_proc_records = 0;
	function<void(FastxRecord_t*, size_t)> callback = process_record;

	if (fastx_ctx.format == FileFormat::FILE_FORMAT_UNKNOWN)
		throw runtime_error("Unknown file format");

	do
	{
		num_proc_records = dispatch_records(in_buf, in_buf_size, 
											&fastx_ctx, 
											&fastx_record,
											&record_idx,
											callback);
	} while (num_proc_records);
}

int64_t get_nth_value(uint64_t col_idx, uint8_t nuc_idx, uint64_t q)
{
	int64_t pos = 0;

	if (col_idx > fastx_ctx.max_seq_len)
		throw out_of_range(format("Invalid range: col_idx={}, nuc_idx={}", col_idx, nuc_idx));

	if (q == 0)
		return col_stats[col_idx].nuc_stats[nuc_idx].min;

	if (q >= col_stats[col_idx].nuc_stats[nuc_idx].count)
		throw out_of_range(format("Invalid range: quantile={}", q));

	while (q > 0)
	{
		if (col_stats[col_idx].nuc_stats[nuc_idx].base_counts[pos] > q)
			break;

		q -= col_stats[col_idx].nuc_stats[nuc_idx].base_counts[pos];
		pos++;

		while (col_stats[col_idx].nuc_stats[nuc_idx].base_counts[pos] == 0)
			pos++;
	}

	return pos + min_qual;
}

void print_headers(Nucleotide nuc = Nucleotide::UNDEFINED)
{
	for (const auto& header : COMMON_HEADERS)
	{
		if (out_ver == OutputVersion::V1)
			fprintf(fastx_ctx.out_stream, "\t%s", header.c_str());
		else
		{
			if (nuc == Nucleotide::ALL)
				fprintf(fastx_ctx.out_stream, "\tALL_%s", header.c_str());
			else
				fprintf(fastx_ctx.out_stream, "\t%c_%s", NUC_CHARS[static_cast<int>(nuc)], header.c_str());
		}
	}
}

void print_nuc_stats(uint64_t col_idx, uint8_t nuc_idx)
{
	int64_t q1, q3, iqr;
	int64_t left_wisker, right_wisker;
	uint64_t i = col_idx;

	q1 = get_nth_value(i, nuc_idx, col_stats[i].nuc_stats[nuc_idx].count / 4);
	q3 = get_nth_value(i, nuc_idx, col_stats[i].nuc_stats[nuc_idx].count * 3 / 4);
	iqr = q3 - q1;

	left_wisker = q1 - iqr * 3 / 2;
	right_wisker = q3 + iqr * 3 / 2;

	if (left_wisker < col_stats[i].nuc_stats[nuc_idx].min)
		left_wisker = col_stats[i].nuc_stats[nuc_idx].min;

	if (right_wisker > col_stats[i].nuc_stats[nuc_idx].max)
		right_wisker = col_stats[i].nuc_stats[nuc_idx].max;

	if (out_ver == OutputVersion::V1)
		fprintf(fastx_ctx.out_stream, "%" PRIu64 "\t", i + 1);

	fprintf(fastx_ctx.out_stream, "%" PRIu64 "\t%d\t%d\t%" PRIu64 "\t",
		col_stats[i].nuc_stats[nuc_idx].count,
		col_stats[i].nuc_stats[nuc_idx].min,
		col_stats[i].nuc_stats[nuc_idx].max,
		col_stats[i].nuc_stats[nuc_idx].sum);

	fprintf(fastx_ctx.out_stream, "%3.2f\t%" PRId64 "\t%" PRId64 "\t%" PRId64 "\t",
		(static_cast<double>(col_stats[i].nuc_stats[nuc_idx].sum)) / (static_cast<double>(col_stats[i].nuc_stats[nuc_idx].count)),
		q1,
		get_nth_value(i, nuc_idx, col_stats[i].nuc_stats[nuc_idx].count / 2),
		q3);

	fprintf(fastx_ctx.out_stream, "%" PRId64 "\t%" PRId64 "\t%" PRId64, iqr, left_wisker, right_wisker);
}

void print_v1_stats()
{
	fprintf(fastx_ctx.out_stream, "column");
	print_headers();
	fprintf(fastx_ctx.out_stream, "\tA_Count\tC_Count\tG_Count\tT_Count\tN_Count\t");
	fprintf(fastx_ctx.out_stream, "Max_count\n");

	for (size_t i = 0; i < fastx_ctx.max_seq_len; ++i)
	{
		if (col_stats[i].nuc_stats[ALL].count == 0)
			break;

		print_nuc_stats(i, ALL);

		fprintf(fastx_ctx.out_stream, "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t",
			col_stats[i].nuc_stats[static_cast<int>(Nucleotide::A)].count,
			col_stats[i].nuc_stats[static_cast<int>(Nucleotide::C)].count,
			col_stats[i].nuc_stats[static_cast<int>(Nucleotide::G)].count,
			col_stats[i].nuc_stats[static_cast<int>(Nucleotide::T)].count,
			col_stats[i].nuc_stats[static_cast<int>(Nucleotide::N)].count);

		fprintf(fastx_ctx.out_stream, "%" PRIu64 "\n", col_stats[0].nuc_stats[ALL].count);
	}
}

void print_v2_stats()
{
	fprintf(fastx_ctx.out_stream, "cycle\tmax_count");

	for (int i = 0; i < NUC_CHARS.size(); ++i)
		print_headers(static_cast<Nucleotide>(i));

	fprintf(fastx_ctx.out_stream, "\n");

	uint64_t max_count = col_stats[0].nuc_stats[static_cast<int>(Nucleotide::ALL)].count;

	for (size_t i = 0; i < fastx_ctx.max_seq_len; ++i)
	{
		if (col_stats[i].nuc_stats[ALL].count == 0)
			break;

		fprintf(fastx_ctx.out_stream, "%zd\t%zd\t", i + 1, max_count);

		for (int j = 0; j < NUC_CHARS.size(); ++j)
			print_nuc_stats(i, static_cast<uint8_t>(j));

		fprintf(fastx_ctx.out_stream, "\n");
	}
}

void print_stats()
{
	if (out_ver == OutputVersion::V1)
		print_v1_stats();
	else
		print_v2_stats();
}

int main(int argc, char** argv)
{
	try
	{
		ios::sync_with_stdio(false);
		parse_args(argc, argv);
		init_vals();

		alloc_bufs();
		open_files();

		read_records();
		print_stats();

		close_files();
		free_bufs();
	}
	catch (const exception& e)
	{
		cout << e.what() << endl;
		exit(errno);
	}

	return 0;
}