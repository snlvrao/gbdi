// g++ -O0 -lm global_base_table.cpp huffman_encoder.cpp main.cpp -o gbdi

#include "global_base_table.h"

GlobalBaseTable global_base_table[NUM_GLOBAL_BASES];
InputDataCount input_data_count;
uint64_t g_compressed_data_size;
uint32_t input_data[NUM_INPUT_ELEMENTS];

int main(int argc, char *argv[])
{
    uint32_t sum = 0;
    uint32_t buffer[128];
    uint64_t block_count = 0;
    uint64_t word_count = 0;
    uint32_t data_count = 0;
    uint32_t num_bins = pow(2, atoi(argv[2]));
    uint64_t original_data_size_mb = 0;
    uint64_t compressed_data_size_mb = 0;

    float_t compression_ratio = 0.0;
    float_t percent_done = 0.0;

    InputDataBlock input_data_block[MAX_DATA_BLOCKS];
    CompressedDataBlock compressed_data_block[MAX_DATA_BLOCKS];

    uint16_t deltas_array[NUM_WORDS_PER_BLOCK];
    uint32_t outliers_array[NUM_WORDS_PER_BLOCK];

    const char *fileName = argv[1];

    ifstream inFile(fileName, ios::binary);

    // ofstream outFile("600.perlbench_s_5_zero_removed.dump", ios::binary);

    while (inFile.read(reinterpret_cast<char *>(&buffer), sizeof(buffer)))
    {
        if(accumulate(buffer, buffer + 128, sum) == 0)
        {
            continue;
        }
        
        for(int i = 0; i < 128; i++)
        {
            input_data_block[block_count].data_word[word_count] = buffer[i];

            if(data_count < NUM_INPUT_ELEMENTS)
            {
                input_data[data_count] = buffer[i];
                data_count++;
            }

            word_count++;
            input_data_count.num_elements++;

            if (word_count == NUM_WORDS_PER_BLOCK)
            {
                block_count++;
                input_data_count.num_blocks++;
                word_count = 0;
            }
        }

        //outFile.write(reinterpret_cast<char *>(&buffer), sizeof(buffer));
    }

    inFile.close();

    original_data_size_mb = (input_data_count.num_elements * sizeof(uint32_t) / (uint32_t) pow(2, 20));
    printf("\noriginal data size = %u MB\n\n", original_data_size_mb);

    compute_global_base_table(input_data_block, global_base_table, num_bins);

    for (uint64_t block_count = 0; block_count < input_data_count.num_blocks; block_count++)
    {
        compute_deltas(&input_data_block[block_count], global_base_table);

        pack_deltas_outliers(&input_data_block[block_count], global_base_table, deltas_array, outliers_array);

        compress_data(&input_data_block[block_count], global_base_table, deltas_array, outliers_array, &compressed_data_block[block_count]);

        percent_done = ((float_t) (block_count + 1) / (float_t) input_data_count.num_blocks) * 100.0;
        printf("Percentage done = %.2f%\r", percent_done);
    }

    compressed_data_size_mb = (uint32_t)(g_compressed_data_size / (uint32_t) (8 * pow(2, 20)));
    
    printf("\n\n[N = 2^%u] num_elements = %u :: total_inlier_count = %u :: total_outlier_count = %u\n", (uint32_t) log2(num_bins), input_data_count.num_elements, input_data_count.total_inlier_count, input_data_count.total_outlier_count);
    printf("compressed_data_size = %u MB\n", compressed_data_size_mb);
    printf("\ncompression ratio = %f\n\n", (float_t) original_data_size_mb / (float_t) compressed_data_size_mb);

    return 0;
}
