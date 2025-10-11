#pragma once

#include <cstdint>
#include <functional>
#include <vector>

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
	std::string in_name;
	std::string out_name;

	FILE* in_stream = stdin;
	FILE* out_stream = stdout;

	FileFormat format = FileFormat::FILE_FORMAT_UNKNOWN;

	size_t max_seq_len = MAX_SEQUENCE_LENGTH;
	size_t total_read_lines = 0;
	size_t total_read_records = 0;
	size_t total_seq_count = 0;
} FastxContext_t;

#pragma pack(push, 1)
typedef struct FastxRecord_s
{
	std::unique_ptr<char[]> seq_id;
	std::unique_ptr<char[]> seq;
	std::unique_ptr<char[]> desc;
	std::unique_ptr<char[]> qual;

	size_t seq_id_len;
	size_t seq_len;
	size_t desc_len;
	size_t qual_len;
	size_t read_count;

	std::unique_ptr<char[]>* operator[](int idx) {
		switch (idx) {
		case 0: return &seq_id;
		case 1: return &seq;
		case 2: return &desc;
		case 3: return &qual;
		default: return nullptr;
		}
	}
} FastxRecord_t;
#pragma pack(pop)

typedef struct DispatchResult_s
{
	size_t num_read_bytes = 0;
	size_t num_proc_bytes = 0;
	size_t num_rem_bytes = 0;
	size_t num_proc_records = 0;
} DispatchResult_t;

static uint8_t RECORD_MEMBER_COUNTS[static_cast<uint8_t>(FileFormat::_FILE_FORMAT_COUNT_)] = { 0, 2, 4 };

static const std::vector<char> FILE_SIGNATURES = { '\0', '>', '@' };

void open_file(const char* file_name, const char* mode, FILE** stream);
void close_file(FILE* stream);

FileFormat get_file_format(FILE* stream);
uint32_t get_read_count(FastxContext_t* ctx, char* seq_id, size_t len);

size_t fwrite_with_line(const void* buf, size_t elem_size, size_t elem_count, FILE* stream, bool use_crlf = false);
DispatchResult_t dispatch_records(char* buf, size_t buf_size, size_t offset, FastxContext_t* ctx, FastxRecord_t* records, size_t* record_idx, std::function<void(FastxRecord_t*, size_t)>& callback);