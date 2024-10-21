#include <cstdio>
#include <cstring>
#include <chrono>
#include <format>
#include <iostream>
#include <stdexcept>
#include "fastx.hpp"

using namespace std;

void open_file(const char* file_name, const char* mode, FILE** stream)
{
	if (file_name && strlen(file_name) > 0)
	{
#ifdef _MSC_VER
		if (fopen_s(stream, file_name, mode))
#else
		* stream = fopen(file_name, mode);
		if (!*stream)
#endif
			throw runtime_error(format("Failed to open file(): {}", errno, file_name));
	}
}

void close_file(FILE* stream)
{
	if (stream && stream != stdin && stream != stdout)
		fclose(stream);
}

FileFormat get_file_format(FILE* stream)
{
	char sig = static_cast<char>(fgetc(stream));
	fseek(stream, -1, SEEK_CUR);

	auto it = find(FILE_SIGNATURES.begin(), FILE_SIGNATURES.end(), sig);
	char idx = it != FILE_SIGNATURES.end() ? static_cast<char>(distance(FILE_SIGNATURES.begin(), it)) : 0;

	return static_cast<FileFormat>(idx);
}

uint32_t get_read_count(FastxContext_t* ctx, char* seq_id, size_t len)
{
	uint32_t result = 1;

	if (!seq_id)
		throw runtime_error("seq_id is null");

	if (ctx->format == FileFormat::FILE_FORMAT_FASTA)
	{
		char* dash = reinterpret_cast<char*>(memchr(seq_id, '-', len));

		if (dash)
			result = atoi(dash + 1);
	}

	return result ? result : 1;
}

size_t fwrite_with_line(const void* buf, size_t elem_size, size_t elem_count, FILE* stream, bool use_crlf)
{
	size_t num_written_bytes = fwrite(buf, elem_size, elem_count, stream);
	num_written_bytes += fwrite(use_crlf ? "\r\n" : "\n", sizeof(char), use_crlf ? 2 : 1, stream);

	return use_crlf ? num_written_bytes + 2 : num_written_bytes + 1;
}

size_t dispatch_lines(char* buf, size_t buf_size, FILE* stream, function<void(const char*, size_t, size_t)> callback)
{
	size_t line_idx = 0;
#ifdef _MSC_VER
	size_t num_read_bytes = fread_s(buf, buf_size, sizeof(char), buf_size, stream);
#else
	size_t num_read_bytes = fread(buf, sizeof(char), buf_size, stream);
#endif
	size_t num_proc_bytes = 0;
	char* prev_eol = buf;
	char* next_eol = nullptr;

	if (num_read_bytes > 0)
	{
		while (prev_eol < buf + num_read_bytes)
		{
			size_t adv_count = 0;
			next_eol = reinterpret_cast<char*>(memchr(prev_eol, LINE_FEED, num_read_bytes - num_proc_bytes));

			if (!next_eol)
				break;

			if (prev_eol + 2 <= next_eol && *(next_eol - 1) == CARRIAGE_RETN)
				adv_count = 1;

			callback(prev_eol, next_eol - prev_eol - adv_count, line_idx);
			num_proc_bytes += next_eol - prev_eol + 1;
			prev_eol = next_eol + 1;
			++line_idx;
		}
	}

	fseek(stream, -static_cast<long>(num_read_bytes - num_proc_bytes), SEEK_CUR);

	return num_proc_bytes;
}

size_t dispatch_records(char* buf, size_t buf_size, FastxContext_t* ctx, FastxRecord_t* records, size_t* record_idx, function<void(FastxRecord_t*, size_t)>& callback)
{
	size_t member_count = RECORD_MEMBER_COUNTS[static_cast<int8_t>(ctx->format)];
	size_t line_offset = 0;
#ifdef _MSC_VER
	size_t num_read_bytes = fread_s(buf, buf_size, sizeof(char), buf_size, ctx->in_stream);
#else
	size_t num_read_bytes = fread(buf, sizeof(char), buf_size, ctx->in_stream);
#endif
	size_t num_proc_bytes = 0;
	size_t num_proc_records = 0;
	char* prev_eol = buf;
	char* next_eol = nullptr;

	if (num_read_bytes > 0)
	{
		while (prev_eol < buf + num_read_bytes)
		{
			char* next_member_buf = nullptr;
			size_t adv_offset = 0;
			size_t len = 0;

			next_eol = reinterpret_cast<char*>(memchr(prev_eol, LINE_FEED, num_read_bytes - num_proc_bytes));

			if (!next_eol)
				break;

			if (prev_eol + 2 <= next_eol && *(next_eol - 1) == CARRIAGE_RETN)
				adv_offset = 1;

			line_offset = ctx->total_read_lines % member_count;
			next_member_buf = reinterpret_cast<char**>(&records[*record_idx])[line_offset];

			if (next_member_buf)
			{
				len = next_eol - prev_eol - adv_offset;

				if (len >= ctx->max_seq_len)
					throw out_of_range(format("Line length out of range: curr: {}, max: {}", len, ctx->max_seq_len));

				copy(prev_eol, prev_eol + len, next_member_buf);
				next_member_buf[len] = '\0';
				reinterpret_cast<size_t*>(reinterpret_cast<char*>(&records[*record_idx]) + sizeof(char*) * member_count)[line_offset] = len;
			}

			if (line_offset == member_count - 1)
			{
				records[*record_idx].read_count = get_read_count(ctx, records[*record_idx].seq_id, records[*record_idx].seq_id_len);
				ctx->total_seq_count += records[*record_idx].seq_len;
				callback(&records[*record_idx], num_proc_records);
				++num_proc_records;
				++ctx->total_read_records;
			}

			num_proc_bytes += next_eol - prev_eol + 1;
			prev_eol = next_eol + 1;
			++ctx->total_read_lines;
		}
	}

	fseek(ctx->in_stream, -static_cast<long>(num_read_bytes - num_proc_bytes), SEEK_CUR);

	return num_proc_records;
}