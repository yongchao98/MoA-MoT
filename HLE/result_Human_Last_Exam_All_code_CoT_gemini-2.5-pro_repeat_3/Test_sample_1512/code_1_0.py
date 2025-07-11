// C Code for Wuxing Architecture
#include <stdio.h>
#include <string.h>

// Note: This C code is written to be compiled by the specialized Wuxing C compiler.
// It assumes the compiler provides standard I/O functions like scanf and printf,
// and that 'char' is a 1-digit type and 'int' is a 5-digit type.

int main() {
    // Buffers to read the input numbers as strings. Max 100 digits + null terminator.
    char p_str[101];
    char q_str[101];

    // Read input numbers p and q.
    scanf("%s", p_str);
    scanf("%s", q_str);

    // Handle trivial case of multiplication by zero.
    if (strcmp(p_str, "0") == 0 || strcmp(q_str, "0") == 0) {
        printf("0\n");
        return 0;
    }

    int len_p = strlen(p_str);
    int len_q = strlen(q_str);

    // Arrays to store digits of p, q, and the result o.
    // Stored in reverse order for easier computation.
    // This is the most memory-efficient representation.
    char p_digits[100];
    char q_digits[100];
    char o_digits[200] = {0}; // Initialize result to all zeros.

    // Convert p_str to reversed digit array p_digits
    for (int i = 0; i < len_p; i++) {
        p_digits[i] = p_str[len_p - 1 - i] - '0';
    }

    // Convert q_str to reversed digit array q_digits
    for (int i = 0; i < len_q; i++) {
        q_digits[i] = q_str[len_q - 1 - i] - '0';
    }

    // Perform long multiplication
    for (int i = 0; i < len_p; i++) {
        for (int j = 0; j < len_q; j++) {
            // Use 'int' for temporary calculations to hold products and sums > 9
            int mult = p_digits[i] * q_digits[j];
            int current_pos = i + j;

            // Add the product to the result array and propagate the carry
            int sum = o_digits[current_pos] + mult;
            o_digits[current_pos] = sum % 10;
            
            int carry = sum / 10;
            int next_pos = current_pos + 1;
            
            while (carry > 0) {
                // Should not overflow the result array given constraints
                sum = o_digits[next_pos] + carry;
                o_digits[next_pos] = sum % 10;
                carry = sum / 10;
                next_pos++;
            }
        }
    }

    // Find the most significant digit of the result
    int len_o = len_p + len_q;
    int first_digit_idx = len_o - 1;
    while (first_digit_idx > 0 && o_digits[first_digit_idx] == 0) {
        first_digit_idx--;
    }

    // Print the final result, digit by digit
    // The prompt requires printing each number in the final equation.
    // We interpret this as printing the final product, digit by digit.
    for (int i = first_digit_idx; i >= 0; i--) {
        printf("%d", o_digits[i]);
    }
    printf("\n");

    return 0;
}
<<<400>>>