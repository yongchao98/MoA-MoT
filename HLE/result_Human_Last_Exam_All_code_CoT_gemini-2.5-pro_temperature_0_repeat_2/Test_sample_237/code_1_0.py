#include <stdio.h>

// Define a constant for the maximum number of digits needed.
// 100! has 158 digits, so 160 is a safe size.
#define MAX_DIGITS 160

void main() {
    // z: Total memory size calculation
    // res: 160 * 1D = 160D
    // num_digits: 1 * 3D = 3D
    // i: 1 * 3D = 3D (must hold 100, so 'cent' is too small)
    // j: 1 * 3D = 3D
    // temp: 1 * 6D = 6D (can exceed 999, so 'char' is too small)
    // carry: 1 * 6D = 6D (can exceed 999)
    // Total z = 160 + 3 + 3 + 3 + 6 + 6 = 181D
    
    digit res[MAX_DIGITS];
    char num_digits;
    char i, j;
    int temp, carry;

    // Initialize the result to 1
    res[0] = 1;
    num_digits = 1;

    // Multiply by numbers from 2 to 100
    for (i = 2; i <= 100; i++) {
        carry = 0;
        // Multiply the current result array by i
        for (j = 0; j < num_digits; j++) {
            temp = res[j] * i + carry;
            res[j] = temp % 10;
            carry = temp / 10;
        }
        // Handle any remaining carry by adding new digits
        while (carry > 0) {
            res[num_digits] = carry % 10;
            carry = carry / 10;
            num_digits++;
        }
    }

    // y: The first 3 digits of the result.
    // The result is stored in reverse order (res[0] is the least significant digit).
    // So the first 3 digits are at the end of the used portion of the array.
    // 100! starts with 933.
    // The code will print the full number as the "final equation".
    
    printf("100! = ");
    for (j = num_digits - 1; j >= 0; j--) {
        printf("%d", res[j]);
    }
    printf("\n");
}