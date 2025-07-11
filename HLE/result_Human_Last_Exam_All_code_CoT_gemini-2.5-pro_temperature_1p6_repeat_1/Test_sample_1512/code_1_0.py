def generate_c_code():
    c_code = """#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// Defines based on the Wuxing architecture
#define CHUNK_SIZE 5
#define P_Q_LEN 100
#define O_LEN 200
#define P_Q_CHUNKS (P_Q_LEN / CHUNK_SIZE) // 100 / 5 = 20
#define O_CHUNKS (O_LEN / CHUNK_SIZE)   // 200 / 5 = 40
#define BASE 100000L                    // 10^5, as a long constant

// Converts a number string into an array of 5-digit integer chunks.
// The chunks are stored with the least significant at index 0.
void stringToChunks(const char* str, int* chunks, int num_chunks) {
    int len = strlen(str);
    int current_pos = len - 1;
    for (int i = 0; i < num_chunks; ++i) {
        int chunk_val = 0;
        long power_of_10 = 1;
        for (int k = 0; k < CHUNK_SIZE && current_pos >= 0; ++k) {
            chunk_val += (str[current_pos] - '0') * power_of_10;
            power_of_10 *= 10;
            current_pos--;
        }
        chunks[i] = chunk_val;
    }
}

// Prints the result stored in an array of chunks.
void printResultChunks(int* chunks, int num_chunks) {
    int start_idx = num_chunks - 1;
    
    // Find the most significant non-zero chunk to avoid printing leading zeros
    while (start_idx > 0 && chunks[start_idx] == 0) {
        start_idx--;
    }

    // Print the first chunk as a standard integer
    printf("%d", chunks[start_idx]);

    // Print the rest of the chunks, padded with leading zeros to 5 digits
    for (int i = start_idx - 1; i >= 0; i--) {
        printf("%05d", chunks[i]);
    }
}

int main() {
    // Allocate space for input strings (100 digits + null terminator)
    char p_str[P_Q_LEN + 1];
    char q_str[P_Q_LEN + 1];

    // Read the two large numbers from standard input
    scanf("%s", p_str);
    scanf("%s", q_str);
    
    // Allocate arrays for the chunked representation of p and q
    int p_chunks[P_Q_CHUNKS] = {0};
    int q_chunks[P_Q_CHUNKS] = {0};
    
    // Convert the input strings to integer chunks
    stringToChunks(p_str, p_chunks, P_Q_CHUNKS);
    stringToChunks(q_str, q_chunks, P_Q_CHUNKS);

    // Intermediate array for the product, using 'long' (10D) to prevent overflow
    // since int(5D) * int(5D) can result in a 10D number.
    long o_temp_chunks[O_CHUNKS] = {0};
    
    // Perform long multiplication using the chunk arrays
    for (int j = 0; j < P_Q_CHUNKS; ++j) {
        if (q_chunks[j] == 0) continue; // Optimization: skip if multiplier chunk is zero
        for (int i = 0; i < P_Q_CHUNKS; ++i) {
            o_temp_chunks[i + j] += (long)p_chunks[i] * q_chunks[j];
        }
    }

    // Final result array of 5D integers
    int o_chunks[O_CHUNKS] = {0};

    // Normalize the intermediate 'long' array into the final 'int' array
    // by propagating carries across chunks.
    long carry = 0;
    for (int k = 0; k < O_CHUNKS; ++k) {
        long current_val = o_temp_chunks[k] + carry;
        o_chunks[k] = current_val % BASE;
        carry = current_val / BASE;
    }
    
    // Print the final equation
    printf("%s * %s = ", p_str, q_str);
    printResultChunks(o_chunks, O_CHUNKS);
    printf("\\n");
    
    return 0;
}
"""
    print(c_code)

generate_c_code()