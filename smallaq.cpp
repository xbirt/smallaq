/*

MIT License

Copyright (c) 2024 Cristian Dumitrescu

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/


#include <iostream>
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

#undef SILENT_OUTPUT

#ifdef __unix
#define fopen_s(pFile,filename,mode) ((*(pFile))=fopen((filename),(mode)))==NULL
#define sscanf_s sscanf
#endif

/* Function to store an unsigned integer in a variable-size, big-endian format
   The output needs to be a buffer of minimum 8 bytes.

   Ideally, this code should be optimized based on platform endianness. */
void store_uint64_varsize(uint64_t value, uint8_t* output, uint8_t* size) {
    // Determine the size needed to store the value
    if (value <= 0b01111111) { // 127, 0x7F
        // One byte: 0xxxxxxx
        output[0] = (uint8_t)value;
        *size = 1;
    }
    else if (value <= 0b0011111111111111) { // 16383, 0x3FFF
        // Two bytes: 10xxxxxx xxxxxxxx
        output[0] = 0x80 | ((value >> 8) & 0x3F); // Set the first two bits to 10 and the next 6 bits to the value
        output[1] = (uint8_t)(value & 0xFF);      // Set the last 8 bits
        *size = 2;
    }
    else if (value <= 0b00011111111111111111111111111111) { // 536870911, 0x1FFFFFFF
        // Four bytes: 110xxxxx xxxxxxxx xxxxxxxx xxxxxxxx
        output[0] = 0xC0 | ((value >> 24) & 0x1F); // Set the first three bits to 110 and the next 5 bits to the value
        output[1] = (uint8_t)((value >> 16) & 0xFF);
        output[2] = (uint8_t)((value >> 8) & 0xFF);
        output[3] = (uint8_t)(value & 0xFF);
        *size = 4;
    }
    else if (value <= 0b0001111111111111111111111111111111111111111111111111111111111111) { // 2305843009213693951ULL, 0x1FFFFFFFFFFFFFFF
        // Eight bytes: 111xxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
        output[0] = 0xE0 | ((value >> 56) & 0x0F); // Set the first three bits to 111 and the next 4 bits to the value
        output[1] = (uint8_t)((value >> 48) & 0xFF);
        output[2] = (uint8_t)((value >> 40) & 0xFF);
        output[3] = (uint8_t)((value >> 32) & 0xFF);
        output[4] = (uint8_t)((value >> 24) & 0xFF);
        output[5] = (uint8_t)((value >> 16) & 0xFF);
        output[6] = (uint8_t)((value >> 8) & 0xFF);
        output[7] = (uint8_t)(value & 0xFF);
        *size = 8;
    }
    else {
        // Value exceeds the maximum representable range
        *size = 0; // Indicate error by setting size to 0
    }
}


// Function to read an unsigned integer from a variable-size encoded buffer
uint64_t read_uint64_varsize(const uint8_t* input, size_t buffer_size, size_t* read_size) {
    if (buffer_size < 1) {
        // Buffer is too small to contain even a single byte
        *read_size = 0;
        return 0;
    }

    uint8_t first_byte = input[0];

    if ((first_byte & 0x80) == 0x00) {
        // One byte: 0xxxxxxx
        *read_size = 1;
        return first_byte;
    }
    else if ((first_byte & 0xC0) == 0x80) {
        // Two bytes: 10xxxxxx xxxxxxxx
        if (buffer_size < 2) {
            // Buffer is too small to contain two bytes
            *read_size = 0;
            return 0;
        }

        *read_size = 2;
        return ((uint64_t)(input[0] & 0x3F) << 8) | input[1];
    }
    else if ((first_byte & 0xE0) == 0xC0) {
        // Four bytes: 110xxxxx xxxxxxxx xxxxxxxx xxxxxxxx
        if (buffer_size < 4) {
            // Buffer is too small to contain four bytes
            *read_size = 0;
            return 0;
        }

        *read_size = 4;
        return ((uint64_t)(input[0] & 0x1F) << 24) |
            ((uint64_t)input[1] << 16) |
            ((uint64_t)input[2] << 8) |
            input[3];
    }
    else if ((first_byte & 0xF0) == 0xE0) {
        // Eight bytes: 111xxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx
        if (buffer_size < 8) {
            // Buffer is too small to contain eight bytes
            *read_size = 0;
            return 0;
        }

        *read_size = 8;
        return ((uint64_t)(input[0] & 0x0F) << 56) |
            ((uint64_t)input[1] << 48) |
            ((uint64_t)input[2] << 40) |
            ((uint64_t)input[3] << 32) |
            ((uint64_t)input[4] << 24) |
            ((uint64_t)input[5] << 16) |
            ((uint64_t)input[6] << 8) |
            input[7];
    }
    else {
        // Invalid encoding
        *read_size = 0;
        return 0;
    }
}

// Function to store data from a source buffer to a destination buffer with less than 8 bits per entry
// is_nucleotide is needed in order to map uracil to thymine internally
bool store_mapped_data(const uint8_t* src, size_t src_len, uint8_t* dst, size_t dst_buffer_size, uint8_t dst_bit_pos, const uint8_t* mapping, size_t mapping_len, size_t* bytes_written, uint8_t* next_bit_pos, bool is_rna) {
    size_t bits_per_entry = (size_t)ceil(log2(mapping_len)); // Calculate bits required per entry
    size_t dst_pos = 0;
    uint8_t bit_pos = dst_bit_pos;

    if (bit_pos >= 8) {
        return false;
    }

    for (size_t i = 0; i < src_len; ++i) {
        // Find the index of the current source value in the mapping array
        uint8_t value_to_look_for = src[i];
        if (is_rna) { // ugly hack to store RNA data using the thymine index for uracil
            if (value_to_look_for == 'U') { 
                value_to_look_for = 'T';
            }
            else if (value_to_look_for == 'u') {
                value_to_look_for = 't';
            }
        }

        uint8_t* mapped_value_position = (uint8_t*)memchr(mapping, value_to_look_for, mapping_len);
        if (!mapped_value_position) {
            // Value not found in the mapping array
            *next_bit_pos = bit_pos;
            return false;
        }
        
        uint8_t mapped_value = (uint8_t)(mapped_value_position - mapping); // there's always less than 8 bits per mapping array, otherwise this makes no sense

        // Write the bits of the index into the destination buffer
        for (size_t bit = 0; bit < bits_per_entry; ++bit) {
            if (dst_pos >= *bytes_written) {
                if (bit_pos == 0) {
                    // Increment the byte count only when a new byte is started
                    (*bytes_written)++;
                }
            }

            if (mapped_value & (1 << (bits_per_entry - 1 - bit))) {
                dst[dst_pos] |= (1 << (7 - bit_pos));
            }
            else {
                dst[dst_pos] &= ~(1 << (7 - bit_pos));
            }

            bit_pos++;
            if (bit_pos == 8) {
                bit_pos = 0;
                dst_pos++;
                if (dst_pos >= dst_buffer_size) {
                    // Check if the destination buffer has enough space
                    *next_bit_pos = bit_pos;
                    return false;
                }
            }
        }
    }

    *bytes_written = dst_pos;
    *next_bit_pos = bit_pos;

    return true;
}

// Function to read mapped data from a source buffer to a destination buffer
bool read_mapped_data(const uint8_t* src, size_t src_len, uint8_t src_bit_pos, uint8_t* dst, size_t dst_len, size_t num_values, const uint8_t* mapping, size_t mapping_len, size_t* bytes_written, size_t* bytes_read, uint8_t* next_bit_pos, bool is_rna) {
    size_t bits_per_entry = (size_t)ceil(log2(mapping_len)); // Calculate bits required per entry
    size_t src_pos = 0;
    uint8_t bit_pos = src_bit_pos;
    size_t values_read = 0;
    *bytes_written = 0;

    if (bit_pos >= 8) {
        return false;
    }

    while (values_read < num_values) {
        // Check if we need to move to the next byte in the source buffer
        if (bit_pos == 8) {
            bit_pos = 0;
            src_pos++;
        }

        // Check if we have enough bytes left in the source buffer
        if (src_pos >= src_len || *bytes_written >= dst_len) {
            return false; // Not enough data in source buffer
        }

        // Read bits corresponding to the bits_per_entry
        uint8_t value = 0;
        for (size_t bit = 0; bit < bits_per_entry; ++bit) {
            // Check if we need to move to the next byte in the source buffer
            if (bit_pos == 8) {
                bit_pos = 0;
                src_pos++;
            }

            if (src[src_pos] & (1 << (7 - bit_pos))) {
                value |= (1 << (bits_per_entry - 1 - bit));
            }

            bit_pos++;
        }

        // Check if the index is within the mapping array bounds
        if (value >= mapping_len) {
            return false; // Invalid index
        }

        // Store the mapped value into the destination buffer
        if (*bytes_written < dst_len) {
            dst[*bytes_written] = mapping[value];
            if (is_rna) { // ugly hack to replace thymine with uracil for RNA
                if (dst[*bytes_written] == 'T') {
                    dst[*bytes_written] = 'U';
                }
                else if (dst[*bytes_written] == 't') {
                    dst[*bytes_written] = 'u';
                }
            }

            (*bytes_written)++;
        }

        values_read++;
    }

    // Check if we need to move to the next byte in the source buffer
    if (bit_pos == 8) {
        bit_pos = 0;
        src_pos++;
    }

    *bytes_read = src_pos; // Calculate total bytes read from source buffer
    *next_bit_pos = bit_pos;

    return true;
}

uint8_t nucleotide_map[] = {'A', 'C', 'G', 'T'};

uint8_t nucleotide_extended_limited_map[] = { 'A', 'C', 'G', 'T', 'N' };

uint8_t nucleotide_extended_map[] = {
    'A', 'C', 'G', 'T',
    'a', 'c', 'g', 't',
    'N',                // ambiguous, any base
    'M',                // ambiguous, adenine (A) or cytosine (C)
    'R',                // ambiguous, A or G
    'Y',                // ambiguous, C or T
    'S',                // ambiguous, G or C
    'W',                // ambiguous, A or T
    'K',                // ambiguous, G or T
    'B',                // ambiguous, C or G or T
    'D',                // ambiguous, A or G or T
    'H',                // ambiguous, A or C or T
    'V',                // ambiguous, A or C or G
};

uint8_t phred33_map[] = {
    '!', '"', '#', '$', '%', '&', '\'', '(', ')', '*',
    '+', ',', '-', '.', '/', '0', '1', '2', '3', '4',
    '5', '6', '7', '8', '9', ':', ';', '<', '=', '>',
    '?', '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
    'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
    'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\',
    ']', '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f',
    'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
    'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
    '{', '|', '}', '~'
};

uint8_t phred64_map[] = {
    '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I',
    'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
    'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']',
    '^', '_', '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
    'h'
};

/*
// Not used: the extended AA map needs the same amount of bits
uint8_t amino_acid_map[] = {
    'A', // Alanine
    'C', // Cysteine
    'D', // Aspartic acid
    'E', // Glutamic acid
    'F', // Phenylalanine
    'G', // Glycine
    'H', // Histidine
    'I', // Isoleucine
    'K', // Lysine
    'L', // Leucine
    'M', // Methionine
    'N', // Asparagine
    'P', // Proline
    'Q', // Glutamine
    'R', // Arginine
    'S', // Serine
    'T', // Threonine
    'V', // Valine
    'W', // Tryptophan
    'Y'  // Tyrosine
};
*/

// The extended AA map needs the same amount of bits as the regular AA map
uint8_t amino_acid_extended_map[] = {
    'A', // Alanine
    'C', // Cysteine
    'D', // Aspartic acid
    'E', // Glutamic acid
    'F', // Phenylalanine
    'G', // Glycine
    'H', // Histidine
    'I', // Isoleucine
    'K', // Lysine
    'L', // Leucine
    'M', // Methionine
    'N', // Asparagine
    'P', // Proline
    'Q', // Glutamine
    'R', // Arginine
    'S', // Serine
    'T', // Threonine
    'V', // Valine
    'W', // Tryptophan
    'Y', // Tyrosine
    'B', // Aspartic acid or Asparagine (ambiguous)
    'Z', // Glutamic acid or Glutamine (ambiguous)
    'X', // Any amino acid (unknown or ambiguous)
    'U', // Selenocysteine
    'O', // Pyrrolysine
    'J', // Leucine or Isoleucine (ambiguous, less common)
    '*'  // Stop codon
};

#define HEADER_IDENTIFIER 0x4653 /* Gets stored as 'SF'. */
#define HEADER_VERSION 0x3130 /* Gets stored as '01'. */

#pragma pack (push, 1)

typedef struct {
    bool is_nucleotide : 1;             // Indicates if the data is a nucleotide sequence
    bool is_extended : 1;               // Indicates if the data is using the extended nucleotide set
    bool is_extended_limited : 1;       // Indicates if the data is using the extended limited nucleotide set (ACGTN)
    bool is_rna : 1;                    // Indicates if the data is RNA (same encoding as DNA, but U instead of T)
    bool has_quality_data : 1;          // Indicates if the data has associated quality scores
    bool is_phred33_quality_data : 1;   // Indicates if the quality data is phred33 (alternatively phred64)
    bool is_unknown_quality_data : 1;   // Indicates if the quality data is unknown (question marks all over)
} sequence_info;

typedef struct {
    uint16_t identifier;
    uint16_t version;
    uint32_t crc32;
    size_t size;
    sequence_info seq_info;
} file_header;

#pragma pack (pop)

bool determine_file_format(const char* filename, sequence_info* info, size_t max_records_to_read) {
    FILE* file;
    int err = fopen_s(&file, filename, "r");
    if (err != 0) {
        return false;
    }

    info->is_nucleotide = true;
    info->is_extended = false;
    info->is_extended_limited = false;
    info->is_rna = false;
    info->has_quality_data = false;
    info->is_phred33_quality_data = false;
    info->is_unknown_quality_data = true;

    char line[1024];
    size_t sequence_records_checked = 0;

    // Determine if the file is FASTA or FASTQ
    while (fgets(line, sizeof(line), file)) {
        if (line[0] == '@') {
            info->has_quality_data = true;
            break;
        }
        else if (line[0] == '>') {
            info->has_quality_data = false;
            break;
        }
        else {
            // no idea what format that might be.
            return false;
        }
    }

    fseek(file, 0, SEEK_SET);

    bool quality_check_done = !info->has_quality_data;
    bool nucleotide_check_done = false;
    while ((!quality_check_done || !nucleotide_check_done) && fgets(line, sizeof(line), file) && sequence_records_checked < max_records_to_read) {
        if (line[0] == '@' || line[0] == '>') {
            sequence_records_checked++;
        }

        if (!nucleotide_check_done) {
            while (!nucleotide_check_done && fgets(line, sizeof(line), file)) {
                if (line[0] == '+') {
                    // time to read quality data
                    break;
                }

                if (line[0] == '@' || line[0] == '>') {
                    sequence_records_checked++;
                    break;
                }

                size_t line_len = strlen(line);
                for (size_t i = 0; i < line_len; ++i) {
                    uint8_t value = line[i];
                    if ((value == '\n' || value == '\r') && i >= line_len - 2) {
                        break;
                    }

                    if (value == 'U' || value == 'u') {
                        info->is_rna = true; // it's not a problem if we enable this for amino acids
                    }
                    
                    bool is_amino_acid = memchr(amino_acid_extended_map, value, sizeof(amino_acid_extended_map)) != NULL;
                    if (value == 'U') {
                        value = 'T';
                    }

                    bool is_nucleotide = memchr(nucleotide_map, value, sizeof(nucleotide_map)) != NULL;

                    bool is_extended_limited_nucleotide = memchr(nucleotide_extended_limited_map, value, sizeof(nucleotide_extended_limited_map)) != NULL;
                    
                    if (value == 'u') {
                        value = 't';
                    }
                    bool is_extended_nucleotide = memchr(nucleotide_extended_map, value, sizeof(nucleotide_extended_map)) != NULL;
                    
                    if (!is_nucleotide && !info->is_extended) {
                        if (!info->is_extended_limited && is_extended_limited_nucleotide) {
                            info->is_extended_limited = true;
                        }

                        if (is_extended_nucleotide && !is_extended_limited_nucleotide) {
                            info->is_extended = true;
                        }
                    }

                    if (is_amino_acid && !is_nucleotide && !is_extended_nucleotide) {
                        info->is_nucleotide = false;
                        info->is_extended = false; // set this correctly
                        info->is_extended_limited = false; // set this correctly
                        info->is_rna = false; // set this correctly
                        nucleotide_check_done = true;
                        break;
                    }

                    if (is_extended_nucleotide && !is_amino_acid) {
                        nucleotide_check_done = true;
                        break;
                    }

                    if (!is_amino_acid && !is_nucleotide && !is_extended_nucleotide) {
                        // don't know what this input is
                        fclose(file);
                        return false;
                    }
                }
            }
        }

        if (line[0] == '+') {
            if (info->has_quality_data) {
                if (!quality_check_done) {
                    while (!quality_check_done && fgets(line, sizeof(line), file)) {
                        if (line[0] == '@' || line[0] == '>') {
                            sequence_records_checked++;
                            break;
                        }

                        size_t line_len = strlen(line);
                        for (size_t i = 0; i < line_len; ++i) {
                            uint8_t value = line[i];
                            if ((value == '\n' || value == '\r') && i >= line_len - 2) {
                                break;
                            }

                            bool is_phred33 = memchr(phred33_map, value, sizeof(phred33_map)) != NULL;
                            bool is_phred64 = memchr(phred64_map, value, sizeof(phred64_map)) != NULL;
                            bool is_unknown = value == '?'; // part of phred33, but some files only have ?
                            if (!is_unknown && info->is_unknown_quality_data) {
                                info->is_unknown_quality_data = false;
                            }

                            if (is_phred33 && !is_phred64 && !is_unknown) {
                                info->is_phred33_quality_data = true;
                                info->is_unknown_quality_data = false;
                                quality_check_done = true;
                                break;
                            }

                            if (!is_phred33 && !is_phred64) {
                                // don't know what this input is
                                fclose(file);
                                return false;
                            }
                        }
                    }
                }
            }
            else {
                // Quality data in fasta file?!
                fclose(file);
                return false;
            }
        }
    }

    fclose(file);
    return true;
}

#define BUFFER_SIZE_INCREMENT 16384

// Encodes an input file with the specified sequence info and saves it into the output file
bool encode_file(const char* input_filename, const char* output_filename, sequence_info* info) {
    FILE* input_file;
    int err = fopen_s(&input_file, input_filename, "r");
    if (err != 0) {
        perror("Error opening input file");
        return false;
    }

    FILE* output_file;
    err = fopen_s(&output_file, output_filename, "wb");
    if (err != 0) {
        perror("Error opening output file");
        fclose(input_file);
        return false;
    }

    char line[1024];
    uint8_t* buffer = (uint8_t*)malloc(BUFFER_SIZE_INCREMENT);
    if (!buffer) {
        perror("Error allocating buffer");
        fclose(input_file);
        fclose(output_file);
        return false;
    }

    file_header header;
    header.identifier = HEADER_IDENTIFIER;
    header.version = HEADER_VERSION;
    header.crc32 = 0;
    header.size = 0;
    header.seq_info = *info;
    fwrite(&header, sizeof(file_header), 1, output_file);

    size_t buffer_size = BUFFER_SIZE_INCREMENT;
    size_t buffer_pos = 0;
    size_t bytes_written = 0;
    size_t values_written = 0;
    uint8_t next_bit_pos = 0;

    const uint8_t* mapping;
    size_t mapping_len;
    if (info->is_nucleotide) {
        if (info->is_extended) {
            mapping = nucleotide_extended_map;
            mapping_len = sizeof(nucleotide_extended_map);
        }
        else if (info->is_extended_limited) {
            mapping = nucleotide_extended_limited_map;
            mapping_len = sizeof(nucleotide_extended_limited_map);
        }
        else {
            mapping = nucleotide_map;
            mapping_len = sizeof(nucleotide_map);
        }
    }
    else {
        mapping = amino_acid_extended_map;
        mapping_len = sizeof(amino_acid_extended_map);
    }

    size_t bits_per_entry_sequence_data = (size_t)ceil(log2(mapping_len)); // Calculate bits required per sequence entry

    const uint8_t* quality_mapping;
    size_t quality_mapping_len;
    size_t bits_per_entry_quality_data;
    if (info->has_quality_data && !info->is_unknown_quality_data) {
        // Process quality data
        quality_mapping = info->is_phred33_quality_data ? phred33_map : phred64_map;
        quality_mapping_len = info->is_phred33_quality_data ? sizeof(phred33_map) : sizeof(phred64_map);
        bits_per_entry_quality_data = (size_t)ceil(log2(quality_mapping_len)); // Calculate bits required per quality entry
    }
    else {
        quality_mapping = NULL;
        quality_mapping_len = 0;
    }

    bool reading_sequence = false;
    bool reading_quality = false;
    size_t line_num = 0;
    while (fgets(line, sizeof(line), input_file)) {
        line_num++;
        size_t line_len = strcspn(line, "\r\n");
        line[line_len] = '\0'; // Remove newline characters

        if (line[0] == '@' || line[0] == '>') {
            // If previously reading sequence or quality, finalize the current buffer
            if (reading_sequence || reading_quality) {
                uint8_t size;
                uint8_t size_buffer[8];
                store_uint64_varsize(values_written, size_buffer, &size);
                fwrite(size_buffer, 1, size, output_file);
                fwrite(buffer, 1, buffer_pos + (next_bit_pos > 0 ? 1 : 0), output_file);
                buffer_pos = 0;
                values_written = 0;
                next_bit_pos = 0;
                reading_sequence = false;
                reading_quality = false;
            }

            // Start processing sequence data
            reading_sequence = true;

            // Write header without '@' or '>' and with null terminator
            fwrite(line + 1, 1, line_len - 1, output_file);
            fputc('\0', output_file);
            continue;
        }

        if (info->has_quality_data && line[0] == '+') {
            // Start reading quality data
            reading_quality = true;
            reading_sequence = false;
            continue;
        }

        if (reading_sequence) {
            if (buffer_pos + (size_t)ceil(line_len * bits_per_entry_sequence_data / 8) >= buffer_size) {
                buffer_size += BUFFER_SIZE_INCREMENT;
                uint8_t *tmp_buffer = (uint8_t*)realloc(buffer, buffer_size);
                if (!tmp_buffer) {
                    perror("Error reallocating buffer");
                    free(buffer);
                    fclose(input_file);
                    fclose(output_file);
                    return false;
                }

                buffer = tmp_buffer;
            }

            if (!store_mapped_data((const uint8_t*)line, line_len, buffer + buffer_pos, buffer_size - buffer_pos, next_bit_pos, mapping, mapping_len, &bytes_written, &next_bit_pos, info->is_rna)) {
                fprintf(stderr, "Error processing sequence data at line %zu, position %zu, offending character '%c'\n", line_num, bytes_written, line[bytes_written]);
                free(buffer);
                fclose(input_file);
                fclose(output_file);
                return false;
            }

            values_written += line_len;
            buffer_pos += bytes_written;
        }

        if (reading_quality && !info->is_unknown_quality_data) {
            // Process quality data
            if (buffer_pos + (size_t)ceil(line_len * bits_per_entry_quality_data / 8) >= buffer_size) {
                buffer_size += BUFFER_SIZE_INCREMENT;
                uint8_t* tmp_buffer = (uint8_t*)realloc(buffer, buffer_size);
                if (!tmp_buffer) {
                    perror("Error reallocating buffer");
                    fclose(input_file);
                    fclose(output_file);
                    return false;
                }

                buffer = tmp_buffer;
            }

            if (!store_mapped_data((const uint8_t*)line, line_len, buffer + buffer_pos, buffer_size - buffer_pos, next_bit_pos, quality_mapping, quality_mapping_len, &bytes_written, &next_bit_pos, false)) {
                fprintf(stderr, "Error processing quality data at line %zu, position %zu, offending character '%c'\n", line_num, bytes_written, line[bytes_written]);
                free(buffer);
                fclose(input_file);
                fclose(output_file);
                return false;
            }

            values_written += line_len;
            buffer_pos += bytes_written;
        }
        else {
            // The quality data only consists of '?'. We don't need to save it.
        }
    }

    // Save the buffer to the output file if there is any remaining data
    if (reading_sequence || reading_quality) {
        uint8_t size;
        uint8_t size_buffer[8];
        store_uint64_varsize(values_written, size_buffer, &size);
        fwrite(size_buffer, 1, size, output_file);
        fwrite(buffer, 1, buffer_pos + (next_bit_pos > 0 ? 1 : 0), output_file);
    }

    free(buffer);
    fclose(input_file);

    size_t file_size = ftell(output_file);
    fseek(output_file, offsetof(file_header, size), SEEK_SET);
    fwrite(&file_size, sizeof(size_t), 1, output_file);
    fclose(output_file);
    return true;
}

bool decode_file(const char* input_filename, size_t line_size, bool fastq_reprint_header) {
    FILE* input_file;
    int err = fopen_s(&input_file, input_filename, "rb");
    if (err != 0) {
        perror("Error opening input file");
        return false;
    }

    // Read the file header
    file_header header;
    if (fread(&header, sizeof(file_header), 1, input_file) != 1) {
        perror("Error reading file header");
        fclose(input_file);
        return false;
    }

    // Validate the header
    if (header.identifier != HEADER_IDENTIFIER) {
        fprintf(stderr, "Invalid file identifier\n");
        fclose(input_file);
        return false;
    }

    // Read and validate version
    char version_str[3];
    version_str[0] = header.version & 0xFF;
    version_str[1] = (header.version >> 8) & 0xFF;
    version_str[2] = '\0';

    int version;
    if (sscanf_s(version_str, "%d", &version) != 1 || version > 1) {
        fprintf(stderr, "Invalid or unsupported file version\n");
        fclose(input_file);
        return false;
    }

    fseek(input_file, 0, SEEK_END);
    size_t file_size = ftell(input_file);
    fseek(input_file, sizeof(file_header), SEEK_SET);

    if (file_size != header.size) {
        fprintf(stderr, "File size does not match header size\n");
        fclose(input_file);
        return false;
    }

    // Determine the mapping arrays
    const uint8_t* seq_mapping;
    size_t seq_mapping_len;
    if (header.seq_info.is_nucleotide) {
        if (header.seq_info.is_extended_limited) {
            seq_mapping = nucleotide_extended_limited_map;
            seq_mapping_len = sizeof(nucleotide_extended_limited_map);
        }
        else if (header.seq_info.is_extended) {
            seq_mapping = nucleotide_extended_map;
            seq_mapping_len = sizeof(nucleotide_extended_map);
        }
        else {
            seq_mapping = nucleotide_map;
            seq_mapping_len = sizeof(nucleotide_map);
        }
    }
    else {
        seq_mapping = amino_acid_extended_map;
        seq_mapping_len = sizeof(amino_acid_extended_map);
    }

    const uint8_t* qual_mapping = NULL;
    size_t qual_mapping_len = 0;
    if (header.seq_info.has_quality_data) {
        if (header.seq_info.is_unknown_quality_data) {
            qual_mapping = NULL;
        }
        else if (header.seq_info.is_phred33_quality_data) {
            qual_mapping = phred33_map;
            qual_mapping_len = sizeof(phred33_map);
        }
        else {
            qual_mapping = phred64_map;
            qual_mapping_len = sizeof(phred64_map);
        }
    }

    uint8_t header_line_buffer[1024];
    uint8_t num_records_buffer[8];
    uint8_t* data_buffer = NULL;
    size_t data_buffer_size = 0;
    uint8_t* output_buffer = (uint8_t*)malloc(line_size + 1);

    // Calculate bit sizes outside the loop
    size_t seq_bit_size = (size_t)ceil(log2(seq_mapping_len));
    size_t qual_bit_size = header.seq_info.has_quality_data && !header.seq_info.is_unknown_quality_data ? (size_t)ceil(log2(qual_mapping_len)) : 0;

    while (fread(header_line_buffer, 1, 1, input_file) == 1) {
        // Read the header line
        size_t header_len = 0;
        while (header_line_buffer[header_len] != '\0' && header_len < 1024) {
            fread(&header_line_buffer[++header_len], 1, 1, input_file);
        }

        if (header_len >= 1024) {
            perror("Encountered header line exceeding 1024 characters.");
            fclose(input_file);
            free(data_buffer);
            free(output_buffer);
            return false;
        }

        header_line_buffer[header_len] = '\0';

#ifndef SILENT_OUTPUT
        printf("%c%s\n", header.seq_info.has_quality_data ? '@' : '>', header_line_buffer);
#endif

        // Read the number of records
        if (fread(num_records_buffer, 1, 8, input_file) != 8) {
            perror("Error reading number of records");
            fclose(input_file);
            free(data_buffer);
            free(output_buffer);
            return false;
        }

        size_t num_records_size;
        size_t num_records = (size_t)read_uint64_varsize(num_records_buffer, 8, &num_records_size);
        fseek(input_file, (long)num_records_size - 8, SEEK_CUR);

        // Calculate buffer size and reallocate if necessary
        size_t total_bits = num_records * (seq_bit_size + qual_bit_size);
        size_t data_bytes_needed = (total_bits + 7) / 8;
        if (data_bytes_needed > data_buffer_size) {
            uint8_t * temp_ptr = (uint8_t*)realloc(data_buffer, data_bytes_needed);
            if (!temp_ptr) {
                if (data_buffer != NULL) free(data_buffer);
                perror("Error reallocating data buffer");
                fclose(input_file);
                free(output_buffer);
                return false;
            }

            data_buffer = temp_ptr;
            data_buffer_size = data_bytes_needed;
        }

        // Read sequence and quality data
        if (fread(data_buffer, 1, data_bytes_needed, input_file) != data_bytes_needed) {
            perror("Error reading data");
            fclose(input_file);
            free(data_buffer);
            free(output_buffer);
            return false;
        }

        // Process sequence data
        size_t bytes_read = 0;
        size_t src_byte_pos = 0;
        size_t bytes_written = 0;
        uint8_t next_bit_pos = 0;
        size_t remaining_records = num_records;
        while (remaining_records > 0) {
            size_t chunk_size = remaining_records < line_size ? remaining_records : line_size;
            if (!read_mapped_data(data_buffer + src_byte_pos, data_bytes_needed - src_byte_pos, next_bit_pos, output_buffer, line_size, chunk_size, seq_mapping, seq_mapping_len, &bytes_written, &bytes_read, &next_bit_pos, header.seq_info.is_rna)) {
                fprintf(stderr, "Error reading sequence data\n");
                fclose(input_file);
                free(data_buffer);
                free(output_buffer);
                return false;
            }

            src_byte_pos += bytes_read;
            output_buffer[bytes_written] = '\0';
#ifndef SILENT_OUTPUT
            printf("%s\n", output_buffer);
#endif
            remaining_records -= chunk_size;
        }

        // Process quality data if present
        if (header.seq_info.has_quality_data) {
#ifndef SILENT_OUTPUT
            printf("+%s\n", fastq_reprint_header ? (char*)header_line_buffer : "");
#endif
            if (header.seq_info.is_unknown_quality_data) {
#ifndef SILENT_OUTPUT
                for (size_t i = 0; i < num_records; i++) {
                    putchar('?');
                    if ((i + 1) % line_size == 0) putchar('\n');
                }
                if (num_records % line_size != 0) putchar('\n');
#endif
            }
            else {
                remaining_records = num_records;
                while (remaining_records > 0) {
                    size_t chunk_size = remaining_records < line_size ? remaining_records : line_size;
                    if (!read_mapped_data(data_buffer + src_byte_pos, data_bytes_needed - src_byte_pos, next_bit_pos, output_buffer, line_size, chunk_size, qual_mapping, qual_mapping_len, &bytes_written, &bytes_read, &next_bit_pos, false)) {
                        fprintf(stderr, "Error reading quality data\n");
                        fclose(input_file);
                        free(data_buffer);
                        free(output_buffer);
                        return false;
                    }

                    src_byte_pos += bytes_read;
                    output_buffer[bytes_written] = '\0';
#ifndef SILENT_OUTPUT
                    printf("%s\n", output_buffer);
#endif
                    remaining_records -= chunk_size;
                }
            }
        }
    }

    free(data_buffer);
    free(output_buffer);
    fclose(input_file);
    return true;
}

void print_usage(const char* executable) {
    printf("Usage:\n");
    printf("  %s identify <filename>\n", executable);
    printf("  %s encode <input_filename> <output_filename>\n", executable);
    printf("  %s decode <input_filename>\n", executable);
    printf("  %s help\n", executable);
}

void print_info(sequence_info seq_info) {
    printf("%s file, containing %s%s%s data%s.\n",
        seq_info.has_quality_data ? "FASTQ" : "FASTA",
        seq_info.is_nucleotide ? (seq_info.is_rna ? "RNA " : "DNA ") : "",
        seq_info.is_nucleotide ? "nucleotide data" : "amino acid data",
        seq_info.is_extended ? " (extended set)" : (seq_info.is_extended_limited ? " (extended limited set)" : ""),
        seq_info.has_quality_data ? (seq_info.is_unknown_quality_data ? " with unknown quality data" : (seq_info.is_phred33_quality_data ? " with phred-33 quality data" : " with phred-64 quality data")) : "");
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        fprintf(stderr, "Error: Missing action\n");
        print_usage(argv[0]);
        return 1;
    }

    const char* action = argv[1];

    if (strcmp(action, "identify") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Error: Missing filename for identify\n");
            print_usage(argv[0]);
            return 1;
        }
        const char* filename = argv[2];
        
        sequence_info seq_info;
        if (determine_file_format(filename, &seq_info, 10)) {
            printf("%s looks like a ", filename);
            print_info(seq_info);
        }
        else {
            printf("Unable to determine the file format.\n");
            return 2;
        }
    }
    else if (strcmp(action, "encode") == 0) {
        if (argc != 4) {
            fprintf(stderr, "Error: Missing input or output filename for encode\n");
            print_usage(argv[0]);
            return 1;
        }
        const char* input_filename = argv[2];
        const char* output_filename = argv[3];
        printf("Will encode from %s to %s\n", input_filename, output_filename);

        sequence_info seq_info;
        if (determine_file_format(input_filename, &seq_info, 10)) {
            printf("The input file will be treated as a ");
            print_info(seq_info);
            if (!encode_file(input_filename, output_filename, &seq_info)) {
                fprintf(stderr, "An error was encountered while encoding.\n");
                return 3;
            }
        }
        else {
            fprintf(stderr, "Unable to determine the file format.\n");
            return 2;
        }
    }
    else if (strcmp(action, "decode") == 0) {
        if (argc != 3) {
            fprintf(stderr, "Error: Missing input filename for decode\n");
            print_usage(argv[0]);
            return 1;
        }
        const char* input_filename = argv[2];
        if (!decode_file(input_filename, 70, true)) {
            fprintf(stderr, "An error was encountered while decoding.\n");
            return 4;
        }
    }
    else if (strcmp(action, "help") == 0) {
        print_usage(argv[0]);

    }
    else {
        fprintf(stderr, "Error: Unknown action '%s'\n", action);
        print_usage(argv[0]);
        return 1;
    }

    return 0;
}

/*
// Old code, used for testing while developing.

#define BYTE_TO_BINARY_PATTERN "%c%c%c%c%c%c%c%c"
#define BYTE_TO_BINARY(byte)  \
  ((byte) & 0x80 ? '1' : '0'), \
  ((byte) & 0x40 ? '1' : '0'), \
  ((byte) & 0x20 ? '1' : '0'), \
  ((byte) & 0x10 ? '1' : '0'), \
  ((byte) & 0x08 ? '1' : '0'), \
  ((byte) & 0x04 ? '1' : '0'), \
  ((byte) & 0x02 ? '1' : '0'), \
  ((byte) & 0x01 ? '1' : '0')

int main()
{
    uint8_t output[8]; // Output buffer
    uint8_t size;      // Variable to store the size

    uint64_t value = 1223428111; // Example value

    store_uint64_varsize(value, output, &size);

    printf("Value stored: %llu\n", value);
    printf("Encoded size: %d bytes\n", size);
    printf("Encoded output: ");
    for (int i = 0; i < size; i++) {
        printf("%02X ", output[i]);
    }
    printf("\n");

    size_t read_size;
    uint64_t read_value = read_uint64_varsize(output, sizeof(output), &read_size);

    printf("Decoded value: %llu\n", read_value);
    printf("Bytes read: %zu\n", read_size);

    printf("------\n");

    uint8_t source_data[] = { 'A', 'T', 'A', 'C', 'G', 'C', 'A', 'G', 'A', 'C', 'G', 'T', 'A', 'C', 'G', 'T' };
    uint8_t encoded_data[sizeof(source_data)];
    memset(encoded_data, 0, sizeof(encoded_data));
    uint8_t initial_dst_bit_pos = 0;
    uint8_t dst_bit_pos = 0;
    size_t bytes_written = 0;

    store_mapped_data(source_data, sizeof(source_data), encoded_data, sizeof(encoded_data), initial_dst_bit_pos, nucleotide_extended_map, sizeof(nucleotide_extended_map), &bytes_written, &dst_bit_pos, false);
    float ratio = 100.0f * (bytes_written * 8 - initial_dst_bit_pos + dst_bit_pos) / (sizeof(source_data) * 8);
    printf("Encoded %zu bytes to %zu bytes (and %d bits). Ratio: %.2f%%.\n", sizeof(source_data), bytes_written, dst_bit_pos, ratio );

    printf("Stored data (hex): ");
    for (int i = 0; i < sizeof(encoded_data); i++) {
        printf("%02X ", encoded_data[i]);
    }

    printf("\n");
    
    printf("Stored data (bin): ");
    for (int i = 0; i < sizeof(encoded_data); i++) {
        printf(BYTE_TO_BINARY_PATTERN, BYTE_TO_BINARY(encoded_data[i]));
    }

    printf("\n");

    uint8_t output_data[sizeof(source_data)];
    memset(output_data, 0, sizeof(output_data));
    size_t values_to_read = sizeof(source_data);
    uint8_t src_bit_pos = 0;
    size_t bytes_read;
    bool is_rna = false;

    read_mapped_data(encoded_data, sizeof(encoded_data), src_bit_pos, output_data, sizeof(output_data), values_to_read, nucleotide_map, sizeof(nucleotide_map), &bytes_written, &bytes_read, &src_bit_pos, is_rna);
    printf("We read %zu bytes, wrote %zu bytes.\n", bytes_read, bytes_written);
    printf("Read data (ascii): ");
    for (size_t i = 0; i < bytes_written; i++) {
        printf("%c", output_data[i]);
    }

    printf("\n");
    printf("------\n");

    sequence_info seq_info;
    const char* input_filename = "SRR28812076.fastq";
    const char* output_filename = "SRR28812076.smallq";

    determine_file_format(input_filename, &seq_info, 10);

    printf("%s file, containing %s%s%s data%s.\n", seq_info.has_quality_data ? "FASTQ" : "FASTA",
        seq_info.is_nucleotide ? (seq_info.is_rna ? "RNA " : "DNA ") : "",
        seq_info.is_nucleotide ? "nucleotide" : "amino acid",
        seq_info.is_extended ? " (extended set)" : (seq_info.is_extended_limited ? " (extended limited set)" : ""),
        seq_info.has_quality_data ? (seq_info.is_unknown_quality_data ? " with unknown quality data" : (seq_info.is_phred33_quality_data ? " with phred-33 quality data" : " with phred-64 quality data")) : "");

    //encode_file(input_filename, output_filename, &seq_info);

    decode_file(output_filename, 70, true);

    return 0;
}
*/
