def factorial_wuxing_program():
    """
    This script generates the C code for calculating 100! on the Wuxing computer
    and prints it in the specified z[C] format.
    """

    # The C code for the Wuxing virtual machine, formatted as a Python string.
    # Note the use of \\n to correctly represent the newline character in the C printf function.
    c_code = """#include <stdio.h>

// This program calculates 100! on the Wuxing computer architecture.
// It is optimized to use the minimum amount of variable memory.
//
// Wuxing C Data Types Used:
// digit: 1D (0-9)
// short: 2D (0-99)
// char:  3D (0-999)

void main() {
    // --- Variable Declaration & Memory Calculation ---
    // digit result[160]: 160 * 1D = 160D
    // char i, j, size:   3 * 3D =   9D
    // char product:      1 * 3D =   3D
    // short carry:        1 * 2D =   2D
    // Total Memory: 160 + 9 + 3 + 2 = 174D

    digit result[160] = {0};
    char i, j, size;
    char product;
    short carry;

    // Initialize result to 1
    result[0] = 1;
    size = 1;

    // Calculate 100! by multiplying from 2 to 100
    for (i = 2; i <= 100; i++) {
        carry = 0;
        // Multiply the current number (in result array) by i
        for (j = 0; j < size; j++) {
            // 'product' fits in char (0-999) because its max value is 9*100+99=999
            product = result[j] * i + carry;
            result[j] = product % 10;
            // 'carry' fits in short (0-99) because its max value is 999/10=99
            carry = product / 10;
        }

        // Propagate the final carry to extend the number's size
        while (carry > 0) {
            result[size] = carry % 10;
            carry = carry / 10;
            size++;
        }
    }

    // Print the final result in the format "100! = [number]"
    printf("100! = ");
    // The result is stored in reverse order, so print from the last digit
    for (j = size - 1; j >= 0; j--) {
        printf("%d", result[j]);
    }
    printf("\\n");
}"""

    # The calculated minimum memory size in 'D' (decimal digits)
    memory_size = 174

    # Print the final answer in the specified z[C] format
    print(f"{memory_size}[{c_code}]")

if __name__ == '__main__':
    factorial_wuxing_program()