def generate_factorial_solution():
    """
    This function generates the C program for calculating 100! on the Wuxing computer
    and calculates its memory footprint as per the problem description.
    """

    # According to the problem, 100! has 158 digits.
    # Data type sizes in 'D' (decimal digits):
    # digit: 1D
    # char: 3D
    #
    # Variable memory calculation:
    # digit result[158]: 158 * 1D = 158D
    # char result_size:   1 * 3D =   3D
    # char i:             1 * 3D =   3D
    # char j:             1 * 3D =   3D
    # char carry:         1 * 3D =   3D
    # Total memory (z): 158 + 3 + 3 + 3 + 3 = 170D
    memory_size = 170

    # The C code is formatted for readability and includes comments
    # explaining the logic and memory usage.
    c_program_code = """void main() {
    // Array to store the digits of the result in reverse order.
    // 100! has 158 digits.
    // Memory size: 158 * sizeof(digit) = 158 * 1D = 158D
    digit result[158];

    // Number of digits currently stored in the result array.
    // Max value will be 158, so char (0-999) is sufficient.
    // Memory size: sizeof(char) = 3D
    char result_size;

    // Loop counter for factorial (i from 2 to 100).
    // Reused for printing loop.
    // Max value 158, so char (0-999) is sufficient.
    // Memory size: sizeof(char) = 3D
    char i;

    // Loop counter for iterating through digits of the result.
    // Max value will be result_size-1 < 158, so char is sufficient.
    // Memory size: sizeof(char) = 3D
    char j;

    // Carry value during multiplication.
    // Max carry: (9 * 100 + 99) / 10 = 99. char (0-999) is sufficient.
    // Memory size: sizeof(char) = 3D
    char carry;

    // --- Algorithm Start ---

    // Initialize result to 1
    result[0] = 1;
    result_size = 1;

    // Loop to calculate 100! by repeated multiplication
    // i goes from 2 up to 100
    for (i = 2; i <= 100; i = i + 1) {
        carry = 0;
        // Multiply the current number (in result array) by i
        for (j = 0; j < result_size; j = j + 1) {
            // Re-use 'carry' variable to store the intermediate product
            carry = result[j] * i + carry;
            result[j] = carry % 10;
            carry = carry / 10;
        }

        // If there is a carry left over, append its digits to the result
        while (carry > 0) {
            result[result_size] = carry % 10;
            carry = carry / 10;
            result_size = result_size + 1;
        }
    }

    // --- Output ---

    // Print the final result.
    // The digits are stored in reverse order in the array.
    // The final equation is 100! = [the big number]
    printf("100! = ");

    // Loop backwards through the array to print digits in correct order.
    // This loop construct (i > 0) works correctly with unsigned char 'i'.
    // This outputs each number (digit) in the final result.
    for (i = result_size; i > 0; i = i - 1) {
        printf("%d", result[i - 1]);
    }
    printf("\\n");
}"""

    # Print the final answer in the format z[C]
    print(f"{memory_size}[{c_program_code}]")

if __name__ == "__main__":
    generate_factorial_solution()