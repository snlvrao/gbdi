#include "global_base_table.h"

bool customer_sorter(BinMetaData const &lhs, BinMetaData const &rhs)
{
    return lhs.count > rhs.count;
}

void compute_global_base_table(InputDataBlock *input_data_block, GlobalBaseTable *global_base_table, uint32_t num_bins)
{
    // ofstream outFile("global_base_table.log", ios::out);

    // TODO: Correct bug in range or find equivalent max_element functionality
    uint64_t max = *max_element(input_data, input_data + NUM_INPUT_ELEMENTS);
    uint64_t min = *min_element(input_data, input_data + NUM_INPUT_ELEMENTS);

    // printf("\nmin = %#08x\n", min);
    // printf("max = %#08x\n\n", max);

    uint32_t bin_size = ceil((double)(max + 1) / num_bins);
    uint32_t input_table[num_bins][bin_size];
    BinMetaData input_bin_count[num_bins];

    u_char global_base_pointer_index_array[NUM_GLOBAL_BASES];
    int frequency[NUM_GLOBAL_BASES];

    // printf("bin_size = %u\n\n", bin_size);

    memset(input_table, 0, sizeof(input_table));
    memset(global_base_table, 0, NUM_GLOBAL_BASES * sizeof(GlobalBaseTable));

    for (int index = 0; index < num_bins; index++)
    {
        input_bin_count[index].count = 0;
        input_bin_count[index].index = index;
    }

    for (int block_index = 0; block_index < (NUM_INPUT_ELEMENTS / NUM_WORDS_PER_BLOCK); block_index++)
    {
        for (int index = 0; index < NUM_WORDS_PER_BLOCK; index++)
        {
            uint32_t data = input_data_block[block_index].data_word[index];
            uint32_t bin_index = (uint32_t) (data / bin_size);
            uint32_t offset = data - (bin_index * bin_size);

            input_table[bin_index][offset]++;
            input_bin_count[bin_index].count++;
        }
    }

    sort(input_bin_count, input_bin_count + num_bins, &customer_sorter);

    for (int index = 0; index < NUM_GLOBAL_BASES; index++)
    {
        uint32_t highest_bin_index = input_bin_count[index].index;
        uint32_t offset = distance(input_table[highest_bin_index], max_element(input_table[highest_bin_index], input_table[highest_bin_index] + bin_size));
        global_base_table[index].base_pointer = (uint32_t)(highest_bin_index * bin_size) + offset;
        global_base_table[index].base_pointer_index = index;
        global_base_table[index].frequency = input_bin_count[index].count;

        global_base_pointer_index_array[index] = index;
        frequency[index] = input_bin_count[index].count;

        // outFile << index << " ";
        // outFile << frequency[index] << endl;

        // printf("value = %#08x\n", global_base_table[index].base_pointer);
        // printf("bin_count = %u\n", input_bin_count[index].count);
        // printf("element_count = %u\n\n", input_table[highest_bin_index][offset]);
    }

    // outFile.close();

    HuffmanCodes(global_base_pointer_index_array, frequency, 150);
}

void compute_deltas(InputDataBlock *input_data_block, GlobalBaseTable *global_base_table)
{
    for (int index = 0; index < NUM_WORDS_PER_BLOCK; index++)
    {
        uint32_t data = input_data_block->data_word[index];
        uint32_t smallest_delta = 0XFFFFFFFF;
        int32_t delta;

        for (int idx = 0; idx < NUM_GLOBAL_BASES; idx++)
        {
            delta = (int32_t)(global_base_table[idx].base_pointer - data);

            if (abs(delta) < smallest_delta)
            {
                smallest_delta = abs(delta);
                input_data_block->delta[index] = delta;
                input_data_block->closest_base_pointer[index] = global_base_table[idx].base_pointer;
                input_data_block->closest_base_pointer_index[index] = idx;
            }
        }

        uint32_t closest_base_pointer_index = input_data_block->closest_base_pointer_index[index];
        delta = abs(input_data_block->delta[index]);

        if ((delta > global_base_table[closest_base_pointer_index].max_delta) && (delta < MAX_DELTA_SIGNED))
        {
            global_base_table[closest_base_pointer_index].max_delta = delta;
        }
        else if (delta >= MAX_DELTA_SIGNED)
        {
            global_base_table[closest_base_pointer_index].max_delta = MAX_DELTA_SIGNED;
        }
        // printf("value = %#08x :: closest base pointer = %#08x :: delta = %#08x\n", input_data_block->data_word[index], input_data_block->closest_base_pointer[index], input_data_block->delta[index]);
    }
}

uint16_t uniqueEle(uint32_t arrayA[], uint32_t arrayB[])
{
   int i,j;
   int count = 0;
   for(i = 0; i < NUM_WORDS_PER_BLOCK; i++)
   {
      for(j = 0; j < i; j++)
      {
         if(arrayA[i] == arrayA[j])
         break;
      }
      if(j == i)
      {
        if(count <= MAX_UNIQUE_BASE_POINTERS)
        {
            arrayB[count] = arrayA[i];
        }
        count++;
      }
   }
   return count;
}

void pack_deltas_outliers(InputDataBlock *input_data_block, GlobalBasetable *global_base_table, uint16_t *deltas_array, uint32_t *outliers_array)
{
    input_data_block->inlier_count = 0;
    input_data_block->outlier_count = 0;
    input_data_block->all_equal = 1;

    for (int index = 0; index < NUM_WORDS_PER_BLOCK; index++)
    {
        uint32_t closest_base_pointer_index = input_data_block->closest_base_pointer_index[index];
        int32_t delta = input_data_block->delta[index];

        if((index > 0) && (input_data_block->data_word[index - 1] != input_data_block->data_word[index]))
        {
            input_data_block->all_equal = 0;
        }

        if ((abs(delta) <= MAX_DELTA_SIGNED))
        {
            input_data_block->mask[index] = 1;
            deltas_array[input_data_block->inlier_count] = (int16_t)delta;
            input_data_block->inlier_count++;
            input_data_count.total_inlier_count++;
            // printf("[inlier] value = %#08x :: closest_base_pointer = %#08x :: delta = %d :: max_delta = %u\n", input_data_block->data_word[index], global_base_table[closest_base_pointer_index].base_pointer, delta, global_base_table[closest_base_pointer_index].max_delta);
        }
        else
        {
            input_data_block->mask[index] = 0;
            outliers_array[input_data_block->outlier_count] = input_data_block->data_word[index];
            input_data_block->outlier_count++;
            input_data_count.total_outlier_count++;
            // printf("[outlier] value = %#08x :: closest_base_pointer = %#08x :: delta = %d :: max_delta = %u\n", input_data_block->data_word[index], global_base_table[closest_base_pointer_index].base_pointer, delta, global_base_table[closest_base_pointer_index].max_delta);
        }
    }
}

void compress_data(InputDataBlock *input_data_block, GlobalBasetable *global_base_table, uint16_t *deltas_array, uint32_t *outliers_array, CompressedDataBlock *compress_data_block)
{
    uint32_t compressed_data_size = 0;
    uint32_t bp_array_len = 0;
    uint32_t deltas_array_len = 0, outliers_array_len = 0;
    uint32_t closest_base_pointer_index = 0;

    input_data_block->num_unique_bp = uniqueEle(input_data_block->closest_base_pointer_index, input_data_block->unique_base_pointer_index);

    // printf("num_unique_bp = %u\n", input_data_block->num_unique_bp);

    // All equal
    if(input_data_block->all_equal)
    {
        compress_data_block->F_enc = 0;
        compress_data_block->uncompressed_val = input_data_block->data_word[0];
        g_compressed_data_size += 32 + F_ENC_LEN;
    }
    else
    {
        if(input_data_block->outlier_count == 0)
        {
            compress_data_block->F_enc = 2;
        }
        else
        {
            compress_data_block->F_enc = 1;
        }

        if(input_data_block->num_unique_bp <= MAX_UNIQUE_BASE_POINTERS)
        {
            compress_data_block->bp_enc = 1;

            bp_array_len = NUM_UNIQUE_BASE_POINTER_BITS; // For #unique base pointers

            for(int i = 0; i < input_data_block->num_unique_bp; i++)
            {
                closest_base_pointer_index = input_data_block->unique_base_pointer_index[i];
                if(global_base_table[closest_base_pointer_index].huffman_used == 1)
                {
                    compress_data_block->H_enc = 1;
                    bp_array_len += global_base_table[closest_base_pointer_index].huffman_code_len;
                }
                else
                {
                    bp_array_len += 11;
                }
            }
        }
        else
        {
            compress_data_block->bp_enc = 0;
        }

        for(int index = 0; index < NUM_WORDS_PER_BLOCK; index++)
        {
            if(input_data_block->mask[index] == 1)
            {
                closest_base_pointer_index = input_data_block->closest_base_pointer_index[index];
                deltas_array_len += (uint32_t) log2(abs(input_data_block->delta[index]));

                if(input_data_block->num_unique_bp <= MAX_UNIQUE_BASE_POINTERS)
                {
                    for(int i = 0; i < input_data_block->num_unique_bp; i++)
                    {
                        if(closest_base_pointer_index == input_data_block->unique_base_pointer_index[i])
                        {
                            compress_data_block->ibp_mask[index] = i;
                        }
                    }
                    bp_array_len += 2;
                }
                else
                {
                    if(global_base_table[closest_base_pointer_index].huffman_used == 1)
                    {
                        compress_data_block->H_enc = 1;
                        bp_array_len += global_base_table[closest_base_pointer_index].huffman_code_len;
                    }
                    else
                    {
                        bp_array_len += 11;
                    }
                }
            }
            else
            {
                outliers_array_len += 32;
            }
        }

        compressed_data_size = bp_array_len + deltas_array_len + outliers_array_len + F_ENC_LEN + BP_ENC_LEN + MASK_LEN + H_ENC_LEN;

        g_compressed_data_size += compressed_data_size;
    }

    // printf("inlier_count = %u :: outlier_count = %u\n", input_data_block->inlier_count,input_data_block->outlier_count);
    // printf("bp_array_len = %u :: deltas_array_len = %u :: outliers_array_len = %u\n", bp_array_len, deltas_array_len, outliers_array_len);
    // printf("original data size = %u bits\n", (8 * NUM_WORDS_PER_BLOCK * sizeof(uint32_t)));
    // printf("compressed_data_size = %u bits\n\n", compressed_data_size);
}