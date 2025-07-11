import textwrap

def solve_wuxing_rask():
    # Plan Step 1 & 2: Define the optimized C code as a string.
    # This C code is optimized for speed on the Wuxing architecture by using
    # a chunk-based multiplication algorithm (base 10000) that leverages
    # the 5D `int` and 10D `long` data types effectively.
    c_code = """
    #include <stdio.h>
    #include <string.h>
    #include <stdlib.h>

    // --- Wuxing Architecture Assumptions ---
    // char: 1 decimal digit (0-9)
    // int:  5 decimal digits (0-99999)
    // long: 10 decimal digits (0-9999999999)
    // Standard C library functions like printf, scanf, strlen are available.

    // --- Optimization Strategy ---
    // We process numbers in chunks of 4 digits. This is the base of our number system.
    // A 4-digit number (0-9999) fits in an `int`.
    // The product of two 4-digit chunks (< 10^8) fits in a `long`.
    // This reduces the number of multiplications from 100*100 to 25*25, speeding up the process.

    #define MAX_DIGITS 100
    #define CHUNK_SIZE 4       // Each `int` chunk holds 4 decimal digits.
    #define BASE 10000         // The base for our calculations (10^CHUNK_SIZE).
    #define P_Q_CHUNKS (MAX_DIGITS / CHUNK_SIZE) // 100/4 = 25 chunks for p and q.
    #define RESULT_CHUNKS (P_Q_CHUNKS * 2)       // 200 digits = 50 chunks for the result.

    /**
     * @brief Parses a digit string into an array of integer chunks in reverse order.
     * e.g., "12345678" -> chunks[0]=5678, chunks[1]=1234
     * @param str The input number as a string.
     * @param chunks The output array of integer chunks.
     * @param num_chunks The size of the chunks array.
     */
    void parse_string_to_chunks(const char *str, int *chunks, int num_chunks) {
        int len = strlen(str);
        int chunk_idx = 0;
        
        for(int i = 0; i < num_chunks; ++i) {
            chunks[i] = 0;
        }

        for (int i = len; i > 0 && chunk_idx < num_chunks; i -= CHUNK_SIZE) {
            char temp_str[CHUNK_SIZE + 1];
            int start = i - CHUNK_SIZE;
            if (start < 0) {
                start = 0;
            }
            int sub_len = i - start;
            strncpy(temp_str, str + start, sub_len);
            temp_str[sub_len] = '\\0';
            chunks[chunk_idx++] = atoi(temp_str);
        }
    }

    int main() {
        char p_str[MAX_DIGITS + 2]; // +2 for newline and null terminator
        char q_str[MAX_DIGITS + 2];

        // Read p and q from the input buffer (e.g., memory location 99999)
        // For standard C, this is stdin.
        scanf("%s", p_str);
        scanf("%s", q_str);

        // Allocate memory for chunks. `int` is 5D, `long` is 10D.
        int p_chunks[P_Q_CHUNKS];
        int q_chunks[P_Q_CHUNKS];
        long o_chunks[RESULT_CHUNKS] = {0}; // Initialize result to zero. Use long for sums.

        // Parse input strings into integer chunks
        parse_string_to_chunks(p_str, p_chunks, P_Q_CHUNKS);
        parse_string_to_chunks(q_str, q_chunks, P_Q_CHUNKS);

        // Perform chunk-based multiplication
        for (int i = 0; i < P_Q_CHUNKS; i++) {
            for (int j = 0; j < P_Q_CHUNKS; j++) {
                if (i + j < RESULT_CHUNKS) {
                    o_chunks[i + j] += (long)p_chunks[i] * q_chunks[j];
                }
            }
        }

        // Propagate carries through the result chunks
        for (int i = 0; i < RESULT_CHUNKS - 1; i++) {
            long carry = o_chunks[i] / BASE;
            o_chunks[i + 1] += carry;
            o_chunks[i] %= BASE;
        }

        // Find the most significant non-zero chunk to avoid leading zeros in output
        int first_chunk_idx = RESULT_CHUNKS - 1;
        while (first_chunk_idx > 0 && o_chunks[first_chunk_idx] == 0) {
            first_chunk_idx--;
        }

        // Print the result `o` to the output buffer (e.g., memory location 99999)
        // For standard C, this is stdout.
        printf("%ld", o_chunks[first_chunk_idx]);
        for (int i = first_chunk_idx - 1; i >= 0; i--) {
            // Print intermediate chunks padded with leading zeros to match CHUNK_SIZE
            printf("%04ld", o_chunks[i]);
        }
        printf("\\n");

        return 0;
    }
    """
    
    # Print the C code for the user
    print("--- Optimized C Program for Wuxing Architecture ---")
    print("```c")
    # textwrap.dedent is used to remove leading whitespace from the string
    print(textwrap.dedent(c_code.strip()))
    print("```")
    print("\n" + "="*50 + "\n")

    # Plan Step 3 & 4: Calculate and explain 'm'
    print("--- Minimized Memory Calculation (m) ---")
    print("The minimized total memory use for p, q, and o is calculated by finding the most efficient way to store them, which is digit-by-digit using the 1D 'char' type.")
    
    # Each number in the final equation is represented here.
    mem_p = 100  # p: max 100 digits * 1D/digit
    mem_q = 100  # q: max 100 digits * 1D/digit
    mem_o = 200  # o: max 200 digits (100+100) * 1D/digit
    
    m = mem_p + mem_q + mem_o
    
    print("\nEquation for m:")
    print(f"m = Memory(p) + Memory(q) + Memory(o)")
    print(f"m = {mem_p}D + {mem_q}D + {mem_o}D = {m}D")

    # Final answer in the specified format
    print(f"\n<<<400>>>")


# Execute the function to generate the solution
solve_wuxing_rask()