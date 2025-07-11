import sys

def generate_c_code():
    """
    This function generates the C code for multiplying two large integers
    on the hypothetical Wuxing architecture.
    """
    c_code = r'''#include <stdio.h>
#include <string.h>

/*
 * Wuxing Architecture Specifications:
 * - Numeric System: Decimal
 * - int: 5D (range 0 to 99999)
 * - long: 10D (range 0 to 9999999999)
 * - char: 1D (range 0 to 9)
 *
 * This program multiplies two large integers, p and q, each up to 100 digits.
 *
 * Optimization Strategy:
 * The program uses schoolbook multiplication with a base of 10000. Each number
 * is represented as an array of 4-digit integers. This base is chosen because
 * the product of two 4-digit numbers fits within a 'long', minimizing the number
 * of chunks and operations needed. The result array 'o' is kept as an 'int'
 * array to minimize memory, with carries handled carefully using temporary 'long'
 * variables during calculation.
 */

// Define constants based on the chosen representation
#define MAX_DIGITS 100
#define BASE 10000       // Each array element is a chunk of 4 digits
#define CHUNK_SIZE 4     // log10(BASE)
#define P_ARR_SIZE (MAX_DIGITS / CHUNK_SIZE) // 100/4 = 25
#define Q_ARR_SIZE (MAX_DIGITS / CHUNK_SIZE) // 100/4 = 25
#define O_ARR_SIZE (P_ARR_SIZE + Q_ARR_SIZE) // 25+25 = 50

/**
 * @brief Converts a number string into an array of integers in the chosen base.
 * The numbers are stored in reverse order (least significant chunk at index 0).
 * @param str The input number as a string.
 * @param arr The output integer array.
 * @param arr_len Pointer to store the resulting length of the array.
 */
void stringToBigInt(const char* str, int* arr, int* arr_len) {
    int str_len = strlen(str);
    *arr_len = 0;
    for (int i = str_len; i > 0; i -= CHUNK_SIZE) {
        int chunk_val = 0;
        int power_of_10 = 1;
        int start = i - CHUNK_SIZE;
        if (start < 0) {
            start = 0;
        }
        for (int k = i - 1; k >= start; k--) {
            chunk_val += (str[k] - '0') * power_of_10;
            power_of_10 *= 10;
        }
        arr[(*arr_len)++] = chunk_val;
    }
}

int main() {
    // Buffers to read the input strings. Max 100 digits + newline + null terminator.
    char p_str[MAX_DIGITS + 2];
    char q_str[MAX_DIGITS + 2];

    // Read p and q from standard input
    scanf("%s", p_str);
    scanf("%s", q_str);

    // Arrays to hold the numbers in base 10000. Initialized to zero.
    int p_arr[P_ARR_SIZE] = {0};
    int q_arr[Q_ARR_SIZE] = {0};
    int p_len, q_len;

    // Convert the input strings to our big integer format
    stringToBigInt(p_str, p_arr, &p_len);
    stringToBigInt(q_str, q_arr, &q_len);

    // Result array 'o'. Using 'int' array for memory optimization.
    int o_arr[O_ARR_SIZE] = {0};

    // Schoolbook multiplication with careful carry handling
    for (int i = 0; i < q_len; i++) {
        // A small optimization: skip if the current chunk of q is zero
        if (q_arr[i] == 0) {
            continue;
        }
        long carry = 0;
        // Multiply the current chunk of q with every chunk of p
        for (int j = 0; j < p_len; j++) {
            // The product and sum can exceed int's max value, so use 'long'
            long prod = (long)q_arr[i] * p_arr[j] + o_arr[i + j] + carry;
            o_arr[i + j] = prod % BASE;
            carry = prod / BASE;
        }
        
        // Propagate the final carry from the inner loop
        int k = i + p_len;
        while (carry > 0 && k < O_ARR_SIZE) {
            long sum = o_arr[k] + carry;
            o_arr[k] = sum % BASE;
            carry = sum / BASE;
            k++;
        }
    }

    // Find the actual length of the resulting number by finding the most significant non-zero chunk
    int final_len = p_len + q_len;
    while (final_len > 1 && o_arr[final_len - 1] == 0) {
        final_len--;
    }

    // Print the final equation with all numbers
    printf("p = %s\n", p_str);
    printf("q = %s\n", q_str);
    printf("o = p * q = ");
    
    // Print the most significant chunk of the result
    printf("%d", o_arr[final_len - 1]);
    
    // Print the remaining chunks, padding with leading zeros to maintain 4 digits
    for (int i = final_len - 2; i >= 0; i--) {
        printf("%04d", o_arr[i]);
    }
    printf("\n");

    return 0;
}
'''
    print(c_code)

if __name__ == "__main__":
    generate_c_code()