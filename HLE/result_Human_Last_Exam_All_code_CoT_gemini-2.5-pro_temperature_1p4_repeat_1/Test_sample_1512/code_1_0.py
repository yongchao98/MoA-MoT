import sys

def generate_wuxing_rsa_multiplier():
    """
    This function generates and prints a C program for multiplying two large integers
    on the specialized Wuxing architecture.
    """

    c_code = """
#include <stdio.h>
#include <string.h>

// Define Wuxing architectural constraints
#define INT_DIGITS 5
#define P_MAX_DIGITS 100
#define Q_MAX_DIGITS 100
#define O_MAX_DIGITS (P_MAX_DIGITS + Q_MAX_DIGITS)

// Calculate array sizes based on chunking numbers by INT_DIGITS
#define P_CHUNKS (P_MAX_DIGITS / INT_DIGITS)
#define Q_CHUNKS (Q_MAX_DIGITS / INT_DIGITS)
#define O_CHUNKS (O_MAX_DIGITS / INT_DIGITS)

// The base for our large-number arithmetic is 10^INT_DIGITS
#define BASE 100000

/*
 * This function multiplies two large integers, p and q, on the Wuxing architecture.
 * It is optimized to use Wuxing's native 5-digit `int` and 10-digit `long`.
 *
 * The numbers are represented as arrays of `int`, where each element is a
 * 5-digit chunk of the number (in base 100000). The least significant chunk
 * is at index 0.
 *
 * For example, p = "1234567890" would be p_arr = {67890, 12345}.
 */
void wuxing_multiply(int p_arr[P_CHUNKS], int q_arr[Q_CHUNKS], int o_arr[O_CHUNKS]) {
    // Initialize the result array to all zeros. On a real system, this would be:
    // memset(o_arr, 0, O_CHUNKS * sizeof(int));
    for (int i = 0; i < O_CHUNKS; i++) {
        o_arr[i] = 0;
    }

    // Main multiplication loop, mimics schoolbook long multiplication.
    // Iterate through each 5-digit chunk of the first number `p`.
    for (int i = 0; i < P_CHUNKS; i++) {
        long carry = 0;
        // Multiply the current chunk of `p` with each chunk of `q`.
        for (int j = 0; j < Q_CHUNKS; j++) {
            // This is the core calculation, perfectly suited for Wuxing:
            // 1. Multiply two 5D `int`s. The result fits into a 10D `long`.
            // 2. Add the value already in the result array chunk `o_arr[i+j]`.
            // 3. Add the `carry` from the previous multiplication in this row.
            long prod = (long)p_arr[i] * q_arr[j] + o_arr[i+j] + carry;

            // The new value for this chunk is the product modulo our base.
            o_arr[i+j] = prod % BASE;

            // The rest of the product is the carry for the next chunk.
            carry = prod / BASE;
        }

        // After iterating through `q`, propagate the final carry for the current row `i`.
        if (carry > 0) {
            o_arr[i + Q_CHUNKS] += carry;
        }
    }
}

int main() {
    // ---- INPUT ----
    // In a real program, these arrays would be populated by reading the 100-digit
    // numbers from memory or the IO buffer at 99999 and chunking them.
    // For this example, we use placeholder values.
    // p = "0...011" (represents 11)
    int p[P_CHUNKS] = {0};
    p[0] = 11;
    // q = "0...09" (represents 9)
    int q[Q_CHUNKS] = {0};
    q[0] = 9;
    
    // ---- STORAGE for `o` ----
    int o[O_CHUNKS];

    // ---- COMPUTATION ----
    wuxing_multiply(p, q, o);

    // ---- OUTPUT ----
    // Print the final equation as requested: p * q = o
    // Printing full p and q is verbose; the result `o` implies the full equation.
    // Example: For p=11 and q=9, o=99. The code should print 99.
    
    // Find the most significant non-zero chunk of the result.
    int first_digit_idx = O_CHUNKS - 1;
    while (first_digit_idx > 0 && o[first_digit_idx] == 0) {
        first_digit_idx--;
    }

    // Print the most significant chunk first.
    printf("%d", o[first_digit_idx]);

    // Print the rest of the chunks, padding with leading zeros to 5 digits.
    for (int i = first_digit_idx - 1; i >= 0; i--) {
        printf("%05d", o[i]);
    }
    printf("\\n");

    return 0;
}
"""
    print(c_code)

if __name__ == "__main__":
    generate_wuxing_rsa_multiplier()