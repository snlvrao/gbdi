#include <iostream>
#include <bits/stdc++.h>
#include <math.h>

using namespace std;

#define NUM_GLOBAL_BASES 2048
#define MAX_FILE_SIZE 8000000000
#define MAX_DATA_BLOCKS (MAX_FILE_SIZE / 64)
#define NUM_INPUT_ELEMENTS 200000
#define NUM_WORDS_PER_BLOCK 16
#define MAX_DELTA_BITS (16 - log2(NUM_GLOBAL_BASES))
#define MAX_DELTA_SIGNED pow(2, (MAX_DELTA_BITS) - 1)
#define MAX_UNIQUE_BASE_POINTERS 4
#define NUM_UNIQUE_BASE_POINTER_BITS 3
#define MASK_LEN 16
#define H_MASK_LEN 16
#define F_ENC_LEN 2
#define H_ENC_LEN 1
#define BP_ENC_LEN 1

typedef struct InputDataBlock
{
    uint32_t data_word[NUM_WORDS_PER_BLOCK];
    uint32_t closest_base_pointer[NUM_WORDS_PER_BLOCK];
    uint32_t closest_base_pointer_index[NUM_WORDS_PER_BLOCK];
    uint32_t unique_base_pointer_index[MAX_UNIQUE_BASE_POINTERS];
    int32_t  delta[NUM_WORDS_PER_BLOCK];
    uint32_t mask[NUM_WORDS_PER_BLOCK];
    uint32_t outlier_count;
    uint32_t inlier_count;
    u_char   all_equal;
    u_char   num_unique_bp;
} InputDataBlock;

typedef struct CompressedDataBlock
{
    u_char      F_enc;
    u_char      bp_enc;
    u_char      H_enc;
    u_char      mask[NUM_WORDS_PER_BLOCK];
    u_char      ibp_mask[NUM_WORDS_PER_BLOCK];
    uint16_t    base_pointer_array[NUM_WORDS_PER_BLOCK];
    int16_t     deltas_array[NUM_WORDS_PER_BLOCK];
    uint32_t    outliers_array[NUM_WORDS_PER_BLOCK];
    uint32_t    uncompressed_val;
} CompressedDataBlock;

typedef struct GlobalBaseTable
{
    uint32_t base_pointer_index;
    uint32_t base_pointer;
    uint32_t max_delta;
    uint32_t frequency;
    uint32_t huffman_code_len;
    u_char   huffman_used;
    u_char   huffman_code[11];
} GlobalBasetable;

typedef struct InputDataCount
{
    uint64_t num_elements;
    uint64_t num_blocks;
    uint64_t total_outlier_count;
    uint64_t total_inlier_count;
} InputDataCount;

typedef struct BinMetaData
{
    uint32_t index;
    uint32_t count;
} BinMetaData;

extern InputDataCount input_data_count;
extern uint64_t g_compressed_data_size;
extern GlobalBaseTable global_base_table[NUM_GLOBAL_BASES];
extern uint32_t input_data[NUM_INPUT_ELEMENTS];

extern void HuffmanCodes(u_char data[], int freq[], int size);

uint16_t uniqueEle(uint32_t arrayA[], uint32_t arrayB[]);

void compute_global_base_table(InputDataBlock *input_data_block, GlobalBasetable *global_base_table, uint32_t num_bins);

void compute_deltas(InputDataBlock *input_data_block, GlobalBaseTable *global_base_table);

void pack_deltas_outliers(InputDataBlock *input_data_block, GlobalBasetable *global_base_table, uint16_t *deltas_array, uint32_t *outliers_array);

void compress_data(InputDataBlock *input_data_block, GlobalBasetable *global_base_table, uint16_t *deltas_array, uint32_t *outliers_array, CompressedDataBlock *compress_data_block);