import sys

def solve():
    """
    This function generates the C code to calculate 100! on the WXVM
    and calculates the minimum memory required for its variables.
    """

    # Step 1 & 4: Calculate minimum memory size 'z'
    # The number of digits in 100! is 158.
    # Memory for the result array: digit result[158]; -> 158 * 1D = 158D
    mem_result = 158
    # Memory for loop counter 'i' (2 to 100): char (0-999) is sufficient. -> 3D
    mem_i = 3
    # Memory for 'result_size' (max 158): char is sufficient. -> 3D
    mem_result_size = 3
    # Memory for inner loop counter 'j' (max 157): char is sufficient. -> 3D
    mem_j = 3
    # Memory for 'carry' (max is < 100): short (0-99) is sufficient. -> 2D
    mem_carry = 2
    # Memory for 'product' (max is 9*100+carry < 999): char is sufficient. -> 3D
    mem_product = 3
    # Memory for 'print_idx' (needs to be signed for reverse loop): int is required. -> 4D
    mem_print_idx = 4
    
    # Total memory z
    z = mem_result + mem_i + mem_result_size + mem_j + mem_carry + mem_product + mem_print_idx

    # Step 2, 3, 5: Write the C code
    c_code = """// WXVM C program to calculate 100! with minimal variable memory.
// It assumes the existence of a stdio.h equivalent for printf.
// The data types (digit, char, short, int) are based on the problem description.

void main() {
    // Variable declarations are optimized for minimum memory usage.
    // Total memory for variables is 176D.
    
    // An array to store the 158 digits of 100!.
    digit result[158]; // 158 * 1D = 158D
    
    // Loop counter for factorial multiplication (from 2 to 100).
    char i; // 3D (range 0-999)
    
    // Tracks the current number of digits in the result array.
    char result_size; // 3D (max value 158)
    
    // Inner loop counter for digit-by-digit multiplication.
    char j; // 3D (max value 157)
    
    // Stores the carry-over during multiplication.
    short carry; // 2D (range 0-99)
    
    // Stores the intermediate product of a digit, 'i', and carry.
    char product; // 3D (range 0-999)
    
    // A signed counter is required for the reverse printing loop.
    int print_idx; // 4D (signed type)

    // Initialize result to 1.
    result[0] = 1;
    result_size = 1;

    // Main loop: Multiply the 'result' array by numbers from 2 to 100.
    for (i = 2; i <= 100; i++) {
        carry = 0;
        for (j = 0; j < result_size; j++) {
            product = result[j] * i + carry;
            result[j] = product % 10;
            carry = product / 10;
        }

        // Append the remaining carry digits to the result array.
        while (carry > 0) {
            result[result_size] = carry % 10;
            result_size++;
            carry = carry / 10;
        }
    }

    // Print the final result. The digits are stored in reverse order.
    // The following loop outputs each digit of the final number, 93326...
    for (print_idx = result_size - 1; print_idx >= 0; print_idx--) {
        printf("%d", result[print_idx]);
    }
}"""

    # Step 6: Format the final output as z[C]
    # The instruction "output each number in the final equation" is interpreted as
    # the C code must print the full value of 100!, which it does.
    final_output = f"{z}[{c_code}]"
    
    # Using sys.stdout.write to avoid adding an extra newline
    sys.stdout.write("<<<")
    sys.stdout.write(final_output)
    sys.stdout.write(">>>")

solve()