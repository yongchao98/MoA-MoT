def solve_factorial_program():
    """
    This script generates the C program for the WXVM to calculate 100!
    and formats it according to the problem specification.
    """

    # The total memory size in 'D' (decimal digits) for the C program's variables.
    memory_size = 494

    # The C program code for the Wuxing Virtual Machine (WXVM).
    c_program_code = """/*
 * Program to calculate 100! on the WXVM.
 *
 * Memory usage for variables (Total: 494D):
 * - result[160]: 160 * sizeof(char) = 160 * 3D = 480D
 * - size:        1 * sizeof(char) = 1 * 3D  =   3D
 * - i, j:        2 * sizeof(char) = 2 * 3D  =   6D
 * - product:     1 * sizeof(char) = 1 * 3D  =   3D
 * - carry:       1 * sizeof(short)= 1 * 2D  =   2D
 */
int main() {
    // Array to hold the digits of the factorial in reverse order.
    // 100! has 158 digits. 160 is a safe size.
    // 'char' (3D, 0-999) is used to store single digits (0-9).
    char result[160] = {0}; 
    
    // 'size' tracks the current number of digits in the result.
    // 'i' and 'j' are loop counters. 'char' is sufficient.
    char size = 1;
    char i, j;
    
    // 'product' holds the intermediate multiplication result for a single digit.
    // Max value is ~ 9 * 100 + 99 = 999. 'char' (0-999) is sufficient.
    char product;

    // 'carry' holds the carry-over for the multiplication.
    // Max value is ~ 999 / 10 = 99. 'short' (0-99) is sufficient.
    short carry;

    // Initialize result to 1 (the value of 1!)
    result[0] = 1;

    // Loop from 2 to 100 to calculate 100!
    for (i = 2; i <= 100; i++) {
        carry = 0;
        // Multiply the number in 'result' array by 'i'
        for (j = 0; j < size; j++) {
            product = result[j] * i + carry;
            result[j] = product % 10;
            carry = product / 10;
        }
        
        // If there's a remaining carry, append its digits to the result array.
        while (carry > 0) {
            result[size] = carry % 10;
            carry = carry / 10;
            size++;
        }
    }
    
    // Print the final result. The digits are stored in reverse order,
    // so we iterate from the end of the array to print the number correctly.
    for (j = size - 1; j >= 0; j--) {
        // %d prints a single digit.
        printf("%d", result[j]);
    }

    return 0;
}"""

    # Print the final answer in the required format: z[C]
    print(f"{memory_size}[{c_program_code}]")

solve_factorial_program()