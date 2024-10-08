#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <format>
#include <iostream>
#include "args.hxx"
#include "fastx.hpp"

using namespace std;

/* Argument parser constants */
static const char* PROLOGUE = "";
static const char* EPILOGUE = "";

/* I/O variables */
static char* in_buf;

/* Argument variables */
static size_t in_buf_size = 32768;
static bool keep_n_nuc_seq = false;
static bool rename_seq_id = false;

/* libfastx variables */
static FastxContext_t fastx_ctx;
static FastxRecord_t fastx_record;

/* internal variables */
static size_t total_bytes_written = 0;

void valid_args()
{
	bool result = true;

	result &= in_buf_size > 0;
	result &= fastx_ctx.max_seq_len > 0;

	if (!result)
		throw invalid_argument("Invalid arguments");
}


void parse_args(int argc, char** argv)
{
	args::ArgumentParser parser(PROLOGUE, EPILOGUE);
	args::HelpFlag help(parser, "help", "Display options", { 'h', "help" });

	args::ValueFlag<string> in_arg(parser, "in", "input file name. default is STDIN", { 'i' }, fastx_ctx.in_name);
	args::ValueFlag<string> out_arg(parser, "out", "output file name. default is STDOUT", { 'o' }, fastx_ctx.out_name);
	args::Flag keep_n_nuc_seq_arg(parser, "n", format("keep sequence with unknown (N) nucleotides. default is {}", keep_n_nuc_seq), { 'n' }, keep_n_nuc_seq);
	args::Flag rn_sqid_arg(parser, "r", format("rename sequence id. default is {}", rename_seq_id), { 'r' }, rename_seq_id);

	args::Group io_tuning_group(parser, "I/O Tuning");
	args::ValueFlag<size_t> ibufs_arg(io_tuning_group, "ibufs", format("input buffer size. default is {}", in_buf_size), { "ibufs" }, in_buf_size);
	args::ValueFlag<size_t> mxsl_arg(io_tuning_group, "mxsl", format("maximum sequence length. default is {}", fastx_ctx.max_seq_len), { "mxsl" }, fastx_ctx.max_seq_len);

	try
	{
		parser.ParseCLI(argc, argv);

		args::get(in_arg).copy(fastx_ctx.in_name, MAX_PATH, 0);
		args::get(out_arg).copy(fastx_ctx.out_name, MAX_PATH, 0);
		keep_n_nuc_seq = args::get(keep_n_nuc_seq_arg);
		rename_seq_id = args::get(rn_sqid_arg);

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

void alloc_bufs()
{
	in_buf = new char[in_buf_size];

	memset(&fastx_record, 0, sizeof(FastxRecord_t));
	fastx_record.seq_id = new char[fastx_ctx.max_seq_len];
	fastx_record.seq = new char[fastx_ctx.max_seq_len];
}

void free_bufs()
{
	if (in_buf) delete[] in_buf;

	if (fastx_record.seq_id) delete[] fastx_record.seq_id;
	if (fastx_record.seq) delete[] fastx_record.seq;
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

void process_record(FastxRecord_t* record, size_t record_idx)
{
	if (!keep_n_nuc_seq && strchr(record->seq, 'N') != NULL)
		return;
	
	if (rename_seq_id)
		fprintf(fastx_ctx.out_stream, "%" PRIu64 "\n", total_bytes_written + 1);
	else
		fprintf(fastx_ctx.out_stream, "%s\n", record->seq_id);

	fprintf(fastx_ctx.out_stream, "%s\n", record->seq);
	total_bytes_written += record->read_count;
}

void read_records()
{
	size_t record_idx = 0; // FIXED: Receive only one record at each time.
	size_t num_proc_records = 0;
	function<void(FastxRecord_t*, size_t)> callback = process_record;

	if (fastx_ctx.format != FileFormat::FILE_FORMAT_FASTQ)
		throw runtime_error("Invalid file format");

	do
	{
		num_proc_records = dispatch_records(in_buf, in_buf_size,
			&fastx_ctx,
			&fastx_record,
			&record_idx,
			callback);
	} while (num_proc_records);
}

int main(int argc, char** argv)
{
	try
	{
		ios::sync_with_stdio(false);
		parse_args(argc, argv);

		alloc_bufs();
		open_files();

		read_records();

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