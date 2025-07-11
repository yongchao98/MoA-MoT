import textwrap

def generate_wuxing_factorial_program():
    """
    This function generates the C code for calculating 100! on the WXVM
    and formats the final output as required.
    """
    
    # z: Total memory size in decimal digits (D)
    # result[200] (digit) = 200 * 1D = 200D
    # i, j, carry, product (char) = 4 * 3D = 12D
    # is_leading_zero (digit) = 1 * 1D = 1D
    # Total = 200 + 12 + 1 = 213D
    z = 213

    # C: The C program code, formatted for clarity.
    # The prompt "output each number in the final equation" is interpreted as
    # printing the digits of the final value of 100!.
    C_code = """\
#include <stdio.h>

/*
 * This program calculates 100! on the WXVM architecture.
 * It is designed to use the smallest amount of memory for variables.
 * Total variable memory usage: 213D.
 */
void main() {
    // --- Variable Declaration (Total 213D) ---
    // The result is stored in reverse order, e.g., 123 as {3, 2, 1, ...}
    digit result[200];      // 200D: Array for the multi-digit result.
    char i;                 //   3D: Loop counter for numbers 2 through 100.
    char j;                 //   3D: Loop counter for traversing the result array.
    char carry;             //   3D: Holds the carry-over during multiplication.
    char product;           //   3D: Temporary storage for (digit * i + carry).
    digit is_leading_zero;  //   1D: Boolean flag (0 or 1) for printing logic.

    // --- Initialization ---
    // Initialize all digits of the result array to 0.
    for (j = 0; j < 200; j++) {
        result[j] = 0;
    }
    // Set the initial value to 1 (for 1!).
    result[0] = 1;

    // --- Factorial Calculation ---
    // Loop from i = 2 to 100, multiplying the result by i in each iteration.
    for (i = 2; i <= 100; i++) {
        carry = 0;
        // Perform long multiplication on the digit array.
        for (j = 0; j < 200; j++) {
            product = result[j] * i + carry;
            result[j] = product % 10;
            carry = product / 10;
        }
    }

    // --- Print Final Result ---
    // This section prints each digit of the calculated value of 100!.
    printf("100! = ");
    
    is_leading_zero = 1;
    // Iterate backwards through the array to print in the correct order.
    for (j = 199; j >= 0; j--) {
        // Start printing only after the most significant digit is found.
        if (result[j] != 0) {
            is_leading_zero = 0;
        }
        if (is_leading_zero == 0) {
            // The '%d' format specifier is for the 'digit' data type.
            printf("%d", result[j]);
        }
    }
    printf("\\n");
}
"""
    
    # Per instructions, the final answer must be in the format z[C]
    # The final code must also output each number in the final equation.
    # My C code interprets this as printing each digit of the result of 100!
    final_output = f"<<<{z}[{textwrap.dedent(C_code)}]\n>>>"
    
    # We print z[C], not the full wrapper, as per outer instructions
    # The wrapper is for the final response block.
    print(f"{z}[{textwrap.dedent(C_code)}]")

# Execute the function to print the final answer.
generate_wuxing_factorial_program()