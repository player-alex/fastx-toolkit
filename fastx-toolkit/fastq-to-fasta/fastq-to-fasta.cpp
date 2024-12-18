#include <cinttypes>
#include <cstdarg>
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
static char* out_buf;

/* Argument variables */
static size_t in_buf_size = 32768;
static size_t out_buf_size = 32768;

static bool keep_n_nuc_seq = false;
static bool rename_seq_id = false;

/* libfastx variables */
static FastxContext_t fastx_ctx;
static FastxRecord_t fastx_record;

/* internal variables */
static uint64_t out_buf_pos = 0;
static string rec_buf;
static size_t total_bytes_written = 0;

void valid_args()
{
	bool result = true;

	result &= in_buf_size > 0 && out_buf_size > 0;
	result &= fastx_ctx.max_seq_len > 0;
	result &= in_buf_size >= fastx_ctx.max_seq_len;

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
	args::ValueFlag<size_t> obufs_arg(io_tuning_group, "obufs", format("output buffer size. default is {}", out_buf_size), { "obufs" }, out_buf_size);
	args::ValueFlag<size_t> mxsl_arg(io_tuning_group, "mxsl", format("maximum sequence length. default is {}", fastx_ctx.max_seq_len), { "mxsl" }, fastx_ctx.max_seq_len);

	try
	{
		parser.ParseCLI(argc, argv);

		args::get(in_arg).copy(fastx_ctx.in_name, MAX_PATH, 0);
		args::get(out_arg).copy(fastx_ctx.out_name, MAX_PATH, 0);
		keep_n_nuc_seq = args::get(keep_n_nuc_seq_arg);
		rename_seq_id = args::get(rn_sqid_arg);

		in_buf_size = args::get(ibufs_arg);
		out_buf_size = args::get(obufs_arg);
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
	out_buf = new char[out_buf_size];

	memset(&fastx_record, 0, sizeof(FastxRecord_t));
	fastx_record.seq_id = new char[fastx_ctx.max_seq_len];
	fastx_record.seq = new char[fastx_ctx.max_seq_len];

	rec_buf.reserve(MAX_SEQUENCE_LENGTH * RECORD_MEMBER_COUNTS[static_cast<uint8_t>(FileFormat::FILE_FORMAT_FASTA)] + 2);
}

void free_bufs()
{
	if (in_buf) delete[] in_buf;
	if (out_buf) delete[] out_buf;

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

void flush_out_buf()
{
	fwrite(out_buf, sizeof(char), out_buf_pos, fastx_ctx.out_stream);
	out_buf_pos = 0;
}

void write_rec(const string& rec)
{
	size_t curr_pos = 0;
	size_t remaining_size = rec.size();

	while (remaining_size > 0)
	{
		size_t copy_size = remaining_size <= out_buf_size ? remaining_size : out_buf_size;

		if (copy_size + out_buf_pos >= out_buf_size)
			copy_size = out_buf_size - out_buf_pos;

		copy(rec.begin() + curr_pos, rec.begin() + curr_pos + copy_size, out_buf + out_buf_pos);

		remaining_size = remaining_size >= copy_size ? remaining_size - copy_size : copy_size - remaining_size;
		out_buf_pos += copy_size;
		curr_pos += copy_size;

		if (out_buf_pos == out_buf_size)
			flush_out_buf();
	}
}

void append_line_break(string& str)
{
	str += LINE_FEED;
}

void process_record(FastxRecord_t* record, size_t record_idx)
{
	if (!keep_n_nuc_seq && strchr(record->seq, 'N') != NULL)
		return;

	if (rename_seq_id)
		rec_buf += to_string(total_bytes_written + 1);
	else
	{
		rec_buf += FILE_SIGNATURES[static_cast<uint8_t>(FileFormat::FILE_FORMAT_FASTA)]; // Append FASTA signature
		rec_buf += record->seq_id + 1; // Skip FASTQ signature
		append_line_break(rec_buf);
	}

	rec_buf += record->seq;
	append_line_break(rec_buf);

	write_rec(rec_buf);
	rec_buf.clear();

	total_bytes_written += record->read_count;
}

void read_records()
{
	size_t record_idx = 0; // FIXED: Receive only one record at each time.
	DispatchResult_t result;
	function<void(FastxRecord_t*, size_t)> callback = process_record;

	if (fastx_ctx.format != FileFormat::FILE_FORMAT_FASTQ)
		throw runtime_error("Invalid file format");

	do
	{
		result = dispatch_records(
			in_buf, in_buf_size,
			result.num_rem_bytes,
			&fastx_ctx,
			&fastx_record,
			&record_idx,
			callback);

	} while (result.num_proc_bytes);

	flush_out_buf();
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