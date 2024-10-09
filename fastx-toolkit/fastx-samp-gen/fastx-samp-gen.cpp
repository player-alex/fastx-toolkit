#include <cstdint>
#include <format>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "args.hxx"
#include "fastx.hpp"

using namespace std;

enum class Nucleotide : unsigned char
{
	ALL, A, C, G, T, N,
	_NUCLEOTIDE_COUNT_
};

static const vector<char> NUC_CHARS = { '\0', 'A', 'C', 'G', 'T', 'N' };

/* Argument parser constants */
static const char* PROLOGUE = "";
static const char* EPILOGUE = "";

/* I/O variables */
static FILE* out_fp = stdout;
static char* out_buf;
static uint64_t out_buf_pos = 0;

/* Argument variables */
static FileFormat file_format = FileFormat::FILE_FORMAT_UNKNOWN;
static int num_records = 0;
static bool collapse_record = false;
static bool use_crlf = false;
static string out_name;

static char base_qual_offset = BASE_QUALITY_OFFSET;
static char min_qual = MIN_QUALITY;
static char max_qual = MAX_QUALITY;

static int min_seq_len = 1;
static int max_seq_len = 50;

static size_t out_buf_size = 32768;

/* Random variables */
static random_device rd;
static mt19937 mt(rd());
static uniform_int_distribution<> alphabet_dist('A', 'Z');
static uniform_int_distribution<> nuc_dist(static_cast<int>(Nucleotide::A), static_cast<int>(Nucleotide::N) - 1);
static uniform_int_distribution<> seq_len_dist;
static uniform_int_distribution<> qual_dist;

void valid_args()
{
	bool result = true;

	result &= !(file_format == FileFormat::FILE_FORMAT_FASTQ && collapse_record);
	result &= num_records > 0;
	
	result &= min_qual < max_qual;
	result &= min_qual < max_qual;

	result &= min_seq_len > 0 && max_seq_len > 0 && min_seq_len <= max_seq_len && max_seq_len <= MAX_SEQUENCE_LENGTH;

	if (!result)
		throw invalid_argument("Invalid arguments");
}

void parse_args(int argc, char** argv)
{
	args::ArgumentParser parser(PROLOGUE, EPILOGUE);

	args::MapFlag<string, FileFormat> sf_arg(parser, "sf", "sample format: fasta or fastq", { 's', "sf" }, {
		{ "fasta", FileFormat::FILE_FORMAT_FASTA },
		{ "fastq", FileFormat::FILE_FORMAT_FASTQ }
	}, args::Options::Required);

	args::ValueFlag<int> nr_arg(parser, "nr", "number of records", { "nr" }, num_records, args::Options::Required);
	args::Flag cr_arg(parser, "cr", "collapse record. only supported to fastq", { 'c' }, collapse_record);
	args::Flag use_crlf_arg(parser, "crlf", format("use crlf. default is {}", use_crlf), { "crlf" }, use_crlf);
	args::ValueFlag<string> out_arg(parser, "out", "output file name. default is STDOUT", { 'o' }, out_name);

	args::Group qual_group(parser, "Quality");
	args::ValueFlag<char> bq_arg(qual_group, "bq", format("base quality offset. default is {}", static_cast<int>(base_qual_offset)), { "bq" }, base_qual_offset);
	args::ValueFlag<char> mnq_arg(qual_group, "mnq", format("min quality. default is {}", static_cast<int>(min_qual)), { "mnq" }, min_qual);
	args::ValueFlag<char> mxq_arg(qual_group, "mxq", format("max quality. default is {}", static_cast<int>(max_qual)), { "mxq" }, max_qual);

	args::Group seq_group(parser, "Sequence");
	args::ValueFlag<int> mns_arg(seq_group, "mns", format("min seq length. default is {}", min_seq_len), { "mns" }, min_seq_len);
	args::ValueFlag<int> mxs_arg(seq_group, "mxs", format("max seq length. max length is {}. default is {}", MAX_SEQUENCE_LENGTH, max_seq_len), { "mxs" }, max_seq_len);

	args::Group io_tuning_group(parser, "I/O Tuning");
	args::ValueFlag<size_t> obufs_arg(io_tuning_group, "obufs", format("output buffer size. default is {}", out_buf_size), { "obufs" }, out_buf_size);

	try
	{
		parser.ParseCLI(argc, argv);

		file_format = args::get(sf_arg);
		num_records = args::get(nr_arg);
		collapse_record = args::get(cr_arg);
		use_crlf = args::get(use_crlf_arg);
		out_name = args::get(out_arg);

		base_qual_offset = args::get(bq_arg);
		min_qual = args::get(mnq_arg);
		max_qual = args::get(mxq_arg);

		min_seq_len = args::get(mns_arg);
		max_seq_len = args::get(mxs_arg);

		out_buf_size = args::get(obufs_arg);

		valid_args();
	}
	catch (exception e)
	{
		cout << parser << endl;
		cerr << e.what() << endl;
		exit(errno);
	}
}

void alloc_bufs()
{
	out_buf = new char[out_buf_size];
}

void free_bufs()
{
	if(out_buf) delete[] out_buf;
}

void append_rand_chars(string& str, size_t len, uniform_int_distribution<>& dist)
{
	for (size_t i = 0; i < len; ++i)
		str += static_cast<char>(dist(mt));
}

void gen_seq_id(string& seq_id, int seq_ord, size_t seq_len)
{	
	seq_id += FILE_SIGNATURES[static_cast<int>(file_format)];
	
	if (file_format == FileFormat::FILE_FORMAT_FASTA)
	{
		seq_id += "sequence";
		seq_id += to_string(seq_ord);
	}

	else if (file_format == FileFormat::FILE_FORMAT_FASTQ)
	{
		append_rand_chars(seq_id, 15, alphabet_dist);
		seq_id += ".";
		seq_id += to_string(seq_ord);

		if (collapse_record)
		{
			append_rand_chars(seq_id, 4, alphabet_dist);
			seq_id += '-';
			append_rand_chars(seq_id, 4, alphabet_dist);
		}

		seq_id += " length=";
		seq_id += to_string(seq_len);
	}
}

void gen_seq(string& seq, size_t seq_len)
{
	for (size_t i = 0; i < seq_len; ++i)
		seq += NUC_CHARS[nuc_dist(mt)];
}

void gen_qual(string& qual, size_t qual_len)
{
	append_rand_chars(qual, qual_len, qual_dist);
}

void flush_out_buf()
{
	fwrite(out_buf, sizeof(char), out_buf_pos, out_fp);
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

		if(out_buf_pos == out_buf_size)
			flush_out_buf();
	}
}

void append_line_break(string& str)
{
	const string new_line = use_crlf ? "\r\n" : "\n";
	str += new_line;
}

void gen_recs()
{
	string rec_buf;
	size_t rec_buf_size = 0;
	uint8_t rec_mem_cnt = RECORD_MEMBER_COUNTS[static_cast<uint8_t>(file_format)];

	rec_buf_size = MAX_SEQUENCE_LENGTH * rec_mem_cnt;
	rec_buf_size += static_cast<size_t>(use_crlf ? 2 : 1) * rec_mem_cnt;
	rec_buf.reserve(rec_buf_size);

	seq_len_dist.param(uniform_int_distribution<>::param_type(min_seq_len, max_seq_len));
	qual_dist.param(uniform_int_distribution<>::param_type(base_qual_offset - min_qual, base_qual_offset + max_qual));

	for (int i = 0; i < num_records; ++i)
	{
		int seq_len = seq_len_dist(mt);

		gen_seq_id(rec_buf, i, seq_len);
		append_line_break(rec_buf);

		gen_seq(rec_buf, seq_len);
		append_line_break(rec_buf);

		if (file_format == FileFormat::FILE_FORMAT_FASTQ)
		{
			rec_buf += '+';
			append_line_break(rec_buf);

			gen_qual(rec_buf, seq_len);
			append_line_break(rec_buf);
		}

		write_rec(rec_buf);
		rec_buf.clear();
	}

	flush_out_buf();
}

int main(int argc, char** argv)
{
	try
	{
		ios::sync_with_stdio(false);
		parse_args(argc, argv);

		open_file(out_name.c_str(), "wb", &out_fp);
		alloc_bufs();

		gen_recs();

		free_bufs();
		close_file(out_fp);
	}
	catch (const exception& e)
	{
		cout << e.what() << endl;
		exit(errno);
	}

	return 0;
}