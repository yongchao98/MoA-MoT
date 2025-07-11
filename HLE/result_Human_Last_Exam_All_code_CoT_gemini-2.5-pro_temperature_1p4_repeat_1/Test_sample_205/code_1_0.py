def generate_wuxing_factorial_program():
    """
    This function constructs the Wuxing C program for calculating 100!
    and formats the output as z[C] where z is memory usage and C is the code.
    """

    # Step 1: Calculate the total memory usage 'z' in decimal digits (D).
    # Based on the problem's data types, we select the smallest possible type for each variable.
    mem_res = 160  # digit res[160]: 160 * 1D
    mem_i = 3      # char i (loops 2 to 100, so short (0-99) is not enough): 1 * 3D
    mem_size = 3   # char size (max ~157): 1 * 3D
    mem_j = 3      # char j (max ~159): 1 * 3D
    mem_carry = 3  # char carry (max product/10 is < 999): 1 * 3D
    total_memory = mem_res + mem_i + mem_size + mem_j + mem_carry

    # Step 2: Define the Wuxing C code as a string.
    # The program calculates 100! and prints the result.
    # The instruction "output each number in the final equation" is fulfilled
    # by the C code's printf("100! = ", ...), which prints the full equation.
    wuxing_c_code = '''#include <stdio.h>

void main() {
    // This program calculates 100! on the Wuxing computer.
    // Total variable memory used: 172D.

    // Variable Declarations:
    digit res[160]; // Array to store factorial digits. 160 * 1D = 160D
    char i;         // Loop counter for 2 to 100. 1 * 3D = 3D
    char size;      // Tracks the number of 'res' digits. 1 * 3D = 3D
    char j;         // Loop counter for digits array. 1 * 3D = 3D
    char carry;     // Stores multiplication carry-over. 1 * 3D = 3D

    // Initialization of the result array
    for (j = 0; j < 160; j = j + 1) {
        res[j] = 0;
    }

    // Initialize result to 1
    res[0] = 1;
    size = 1;

    // Main loop to calculate 100!
    for (i = 2; i <= 100; i = i + 1) {
        carry = 0;
        // Multiply each digit of the result by i
        for (j = 0; j < size; j = j + 1) {
            carry = res[j] * i + carry;
            res[j] = carry % 10;
            carry = carry / 10;
        }
        // Handle remaining carry by adding new digits
        while (carry > 0) {
            res[size] = carry % 10;
            size = size + 1;
            carry = carry / 10;
        }
    }

    // Print the final equation: 100! = [number]
    printf("100! = ");
    
    // Use a while loop for printing as 'char' is unsigned (0-999)
    j = size;
    while (j > 0) {
        j = j - 1;
        printf("%d", res[j]);
    }
    printf("\\n");
}'''

    # Step 3: Print the final answer in the specified z[C] format.
    print(f"{total_memory}[{wuxing_c_code}]")

generate_wuxing_factorial_program()