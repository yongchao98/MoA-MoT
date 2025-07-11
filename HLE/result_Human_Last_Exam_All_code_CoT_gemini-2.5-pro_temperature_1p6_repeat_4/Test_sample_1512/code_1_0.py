import sys
import io

# We need to create a C program as a string and then print it.
# The user wants an optimized C program for the Wuxing architecture to multiply two 100-digit integers.
# Then, we calculate the minimized total memory use for storing p, q, and the result o.

c_program_code = """
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*
 * This is an optimized C program for the hypothetical Wuxing architecture.
 *
 * Wuxing architecture specifications:
 * - Numeric System: Decimal
 * - int: 5 decimal digits (range 0-99999)
 * - long: 10 decimal digits
 */

// Constants based on Wuxing architecture and problem size
#define MAX_DIGITS 100
#define CHUNK_SIZE 5
#define BASE 100000L

#define P_CHUNKS (MAX_DIGITS / CHUNK_SIZE) // 100 / 5 = 20
#define Q_CHUNKS (MAX_DIGITS / CHUNK_SIZE) // 100 / 5 = 20
#define O_CHUNKS (P_CHUNKS + Q_CHUNKS)   // 20 + 20 = 40

// Parses a string of digits into an array of 5-digit integer chunks.
// The least significant chunk is at index 0.
void string_to_chunks(const char* str, int* chunks, int num_chunks) {
    char buf[CHUNK_SIZE + 1];
    buf[CHUNK_SIZE] = '\\0';
    int len = strlen(str);

    for (int i = 0; i < num_chunks; ++i) {
        int start = len - (i + 1) * CHUNK_SIZE;
        int num_copy = CHUNK_SIZE;

        if (start < 0) {
            num_copy += start;
            start = 0;
        }

        if (num_copy > 0) {
            strncpy(buf, str + start, num_copy);
            buf[num_copy] = '\\0'; // Null-terminate the copied part
            chunks[i] = atoi(buf);
        } else {
            chunks[i] = 0;
        }
    }
}

// Main program entry point
int main() {
    // Buffers to read the 100-digit numbers as strings.
    // +2 for potential newline and null terminator.
    char p_str[MAX_DIGITS + 2];
    char q_str[MAX_DIGITS + 2];

    // On a real Wuxing system, I/O might be memory-mapped.
    // For this standard C program, we read from stdin.
    // Example: Two 100-digit numbers on separate lines.
    if (!fgets(p_str, sizeof(p_str), stdin) || !fgets(q_str, sizeof(q_str), stdin)) {
        return 1; // Input error
    }
    p_str[strcspn(p_str, "\\r\\n")] = 0;
    q_str[strcspn(q_str, "\\r\\n")] = 0;

    // Allocate arrays for p, q, and o using Wuxing's 'int' type.
    int p[P_CHUNKS] = {0};
    int q[Q_CHUNKS] = {0};
    int o[O_CHUNKS] = {0};

    // Convert input strings to integer chunk arrays.
    string_to_chunks(p_str, p, P_CHUNKS);
    string_to_chunks(q_str, q, Q_CHUNKS);

    // Temporary storage for intermediate products.
    // Must be 'long' (10D) to hold (p_chunk * q_chunk), which can be ~10^10.
    long temp_o[O_CHUNKS + 1] = {0};

    // Optimized multiplication using 5-digit chunks.
    for (int i = 0; i < Q_CHUNKS; i++) {
        for (int j = 0; j < P_CHUNKS; j++) {
            // Cast to 'long' before multiplication to prevent overflow.
            temp_o[i + j] += (long)p[j] * q[i];
        }
    }

    // Propagate carries through the temporary result array.
    for (int i = 0; i < O_CHUNKS; i++) {
        if (temp_o[i] >= BASE) {
            temp_o[i + 1] += temp_o[i] / BASE;
            temp_o[i] %= BASE;
        }
    }

    // Copy the final chunks from 'long' array to 'int' result array 'o'.
    for (int i = 0; i < O_CHUNKS; i++) {
        o[i] = (int)temp_o[i];
    }

    // Print the final result (up to 200 digits).
    // Find the most significant chunk to avoid printing leading zeros.
    int start_idx = O_CHUNKS - 1;
    while (start_idx > 0 && o[start_idx] == 0) {
        start_idx--;
    }

    // Print the first non-zero chunk without padding.
    printf("%d", o[start_idx]);

    // Print remaining chunks, zero-padded to 5 digits.
    for (int i = start_idx - 1; i >= 0; i--) {
        printf("%05d", o[i]);
    }
    printf("\\n");

    return 0;
}
"""

print("--- Optimized C Program for Wuxing Architecture ---")
print(c_program_code)
print("-------------------------------------------------")
print("\n--- Memory Calculation (m) ---")

# Calculate memory usage in decimal digits (D)
# The minimized memory is achieved by using the most compact data types.
# int (5D), long (10D), char(1D). Storing a 100-digit number requires 100D regardless
# of whether it's char[100], int[20], or long[10].

# Memory for p (100 digits)
mem_p = 100 # D

# Memory for q (100 digits)
mem_q = 100 # D

# Memory for o (up to 200 digits)
mem_o = 200 # D

# Total minimized memory m
m = mem_p + mem_q + mem_o

print(f"Minimized memory to store p (100 digits): {mem_p} D")
print(f"Minimized memory to store q (100 digits): {mem_q} D")
print(f"Minimized memory to store o (200 digits): {mem_o} D")
print("\nThe final equation for the minimized total memory use (m) is:")
print(f"m = {mem_p} + {mem_q} + {mem_o} = {m}")
print("--------------------------------")

sys.stdout = io.StringIO() # Suppress further output in the final block
final_answer = m
sys.stdout = sys.__stdout__ # Restore stdout
# The final answer needs to be enclosed in <<<>>>
# <<<400>>>