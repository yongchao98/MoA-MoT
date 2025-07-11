/*
 * This is a C program. It is placed in a python code block
 * per the instructions.
 */
#include <stdio.h>
#include <string.h>

// A utility function to reverse a string in-place
void reverse_string(char *str) {
    if (str == NULL) {
        return;
    }
    int len = strlen(str);
    int i;
    for (i = 0; i < len / 2; i++) {
        char temp = str[i];
        str[i] = str[len - i - 1];
        str[len - i - 1] = temp;
    }
}

int main() {
    // p and q can have up to 100 digits. We need 101 chars for string input.
    // Memory for p: char p_str[101] -> 101 D
    // Memory for q: char q_str[101] -> 101 D
    char p_str[101];
    char q_str[101];

    // Read inputs using the I/O buffer at location 99999
    scanf("%s", p_str);
    scanf("%s", q_str);

    int len_p = strlen(p_str);
    int len_q = strlen(q_str);

    // The result 'o' can have up to 200 digits.
    // A char array of size 201 holds the result digits (0-9) safely.
    // Memory for o: char o_digits[201] -> 201 D
    char o_digits[201] = {0}; // Initialize with integer 0s

    int i, j;
    
    // Reverse strings to make indexing p[i] as the 10^i coefficient
    reverse_string(p_str);
    reverse_string(q_str);

    // Perform multiplication
    for (i = 0; i < len_q; i++) {
        // According to Wuxing, char is 0-9, so use 'int' for calculations
        int carry = 0;
        int digit_q = q_str[i] - '0';
        
        for (j = 0; j < len_p; j++) {
            int digit_p = p_str[j] - '0';
            
            // Perform the multiplication and add to the existing digit, plus carry
            int product = o_digits[i + j] + (digit_p * digit_q) + carry;
            
            // The new digit is the remainder, the new carry is the quotient
            o_digits[i + j] = product % 10;
            carry = product / 10;
        }
        
        // After iterating through p, propagate the remaining carry
        if (carry > 0) {
            o_digits[i + j] += carry;
            // Handle potential cascade if the sum exceeds 9
            int k = i + j;
            while (o_digits[k] >= 10) {
                int temp_carry = o_digits[k] / 10;
                o_digits[k] = o_digits[k] % 10;
                k++;
                o_digits[k] += temp_carry;
            }
        }
    }
    
    // Find the actual length of the result string by trimming leading zeros
    int len_o = len_p + len_q;
    while (len_o > 1 && o_digits[len_o - 1] == 0) {
        len_o--;
    }

    // Restore original order of p and q for printing
    reverse_string(p_str);
    reverse_string(q_str);

    // Print the final equation
    printf("%s * %s = ", p_str, q_str);
    for (i = len_o - 1; i >= 0; i--) {
        printf("%d", o_digits[i]);
    }
    printf("\n");
    
    return 0;
}