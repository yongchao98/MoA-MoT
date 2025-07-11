import textwrap

# This Python script generates the optimized C code for the Wuxing architecture.
# The user can execute this script to get the C source file.

# The C code implements large number multiplication using a chunking method,
# which is optimized for the Wuxing architecture's 5-digit `int` type.
c_program_string = r"""
#include <stdio.h>
#include <string.h>

// Define constants based on Wuxing architecture and problem constraints
#define MAX_DIGITS 100
#define CHUNK_SIZE 5        // Wuxing `int` size is 5 digits
#define BASE 100000L        // We perform calculations in base 10^5

// Calculate array sizes based on chunks
#define P_CHUNKS (MAX_DIGITS / CHUNK_SIZE) // 20 chunks for p
#define Q_CHUNKS (MAX_DIGITS / CHUNK_SIZE) // 20 chunks for q
#define O_CHUNKS (P_CHUNKS + Q_CHUNKS)   // 40 chunks for the result o

/**
 * @brief Parses a string of digits into an array of 5-digit integer chunks.
 * 
 * The chunks are stored with the most significant at index 0. The input string is
 * processed from right to left to fill the chunk array from least to most significant.
 * 
 * @param s The input string of digits.
 * @param chunks The output array of integers (5D `int`s).
 * @param total_chunks The total size of the chunks array.
 */
void parse_to_chunks(const char* s, int* chunks, int total_chunks) {
    int len = strlen(s);
    
    // Initialize all chunks to 0.
    for (int i = 0; i < total_chunks; i++) {
        chunks[i] = 0;
    }

    int current_chunk_idx = total_chunks - 1;
    for (int i = len; i > 0; i -= CHUNK_SIZE) {
        if (current_chunk_idx < 0) break;

        int start_pos = i - CHUNK_SIZE;
        if (start_pos < 0) {
            start_pos = 0;
        }
        int sub_len = i - start_pos;
        
        // A self-contained routine to convert a portion of a string to an integer.
        int val = 0;
        for (int k = 0; k < sub_len; ++k) {
            val = val * 10 + (s[start_pos + k] - '0');
        }
        chunks[current_chunk_idx--] = val;
    }
}

int main() {
    // We assume input is from stdin. The Wuxing environment would map its
    // I/O buffer (memory location 99999) to this standard stream.
    char p_str[MAX_DIGITS + 2]; // +2 for potential newline and null terminator
    char q_str[MAX_DIGITS + 2];

    // Read the two large numbers as strings.
    scanf("%s", p_str);
    scanf("%s", q_str);
    
    // Data structures for multiplication.
    // Wuxing `int` is 5D, perfect for our chunks.
    int p[P_CHUNKS];
    int q[Q_CHUNKS];
    
    // Wuxing `long` is 10D, necessary for intermediate products to avoid overflow.
    // A product of two 5D numbers (e.g., 99999*99999) is ~10^10, which fits in a 10D long.
    long o_long[O_CHUNKS] = {0};

    // Parse input strings into arrays of 5-digit chunks.
    parse_to_chunks(p_str, p, P_CHUNKS);
    parse_to_chunks(q_str, q, Q_CHUNKS);
    
    // Perform multiplication using the chunk arrays (base 100000).
    // This is computationally more efficient than a simple digit-by-digit approach.
    for (int i = P_CHUNKS - 1; i >= 0; i--) {
        for (int j = Q_CHUNKS - 1; j >= 0; j--) {
            // The product of chunk p[i] and q[j] adds to the result at position o_long[i+j+1].
            o_long[i + j + 1] += (long)p[i] * q[j];
        }
    }

    // Propagate carries through the result array.
    for (int i = O_CHUNKS - 1; i > 0; i--) {
        if (o_long[i] >= BASE) {
            o_long[i - 1] += o_long[i] / BASE;
            o_long[i] %= BASE;
        }
    }

    // Print the final product, o.
    int first_chunk_printed = 0;
    for (int i = 0; i < O_CHUNKS; i++) {
        // Find the first non-zero chunk to begin printing.
        if (o_long[i] != 0 || first_chunk_printed) {
            if (first_chunk_printed) {
                // Print subsequent chunks with leading zeros to ensure they are 5 digits long.
                printf("%05d", (int)o_long[i]);
            } else {
                // Print the most significant chunk without modification.
                printf("%d", (int)o_long[i]);
                first_chunk_printed = 1;
            }
        }
    }

    // If no chunks were printed, the result must be 0.
    if (!first_chunk_printed) {
        printf("0");
    }
    printf("\n");

    return 0;
}
"""

print(textwrap.dedent(c_program_string).strip())