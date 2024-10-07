#pragma once

#include <cstdint>
#include <functional>
#include <vector>

using namespace std;

constexpr char LINE_FEED = '\n';
constexpr char CARRIAGE_RETN = '\r';
constexpr size_t MAX_PATH = 255;
constexpr size_t MAX_SEQUENCE_LENGTH = 25000;

constexpr int BASE_QUALITY_OFFSET = 33;
constexpr int MIN_QUALITY = -15;
constexpr int MAX_QUALITY = 93;

enum class FileFormat : uint8_t
{
	FILE_FORMAT_UNKNOWN,
	FILE_FORMAT_FASTA,
	FILE_FORMAT_FASTQ,
	_FILE_FORMAT_COUNT_,
};

typedef struct FastxContext_s
{
	char in_name[MAX_PATH] = { 0 };
	char out_name[MAX_PATH] = { 0 };

	FILE* in_stream = stdin;
	FILE* out_stream = stdout;

	FileFormat format = FileFormat::FILE_FORMAT_UNKNOWN;

	size_t max_seq_len = MAX_SEQUENCE_LENGTH;
	size_t total_read_lines = 0;
	size_t total_read_records = 0;
} FastxContext_t;

typedef struct FastxRecord_s
{
	char* seq_id;
	char* seq;
	char* desc;
	char* qual;

	size_t seq_id_len;
	size_t seq_len;
	size_t desc_len;
	size_t qual_len;
	size_t read_count;
} FastxRecord_t;

static uint8_t RECORD_MEMBER_COUNTS[static_cast<uint8_t>(FileFormat::_FILE_FORMAT_COUNT_)] = { 0, 2, 4 };

const vector<char> FILE_SIGNATURES = { '\0', '>', '@' };

void open_file(const char* file_name, const char* mode, FILE** stream);
void close_file(FILE* stream);

FileFormat get_file_format(FILE* stream);
uint32_t get_read_count(FastxContext_t* ctx, char* seq_id, size_t len);

size_t fwrite_with_line(const void* buf, size_t elem_size, size_t elem_count, FILE* stream);
size_t dispatch_lines(char* buf, size_t buf_size, FILE* stream, function<void(const char*, size_t, size_t)> callback);

size_t dispatch_records(char* buf, size_t buf_size, FastxContext_t* ctx, FastxRecord_t* records, size_t* record_idx, function<void(FastxRecord_t*, size_t)>& callback);