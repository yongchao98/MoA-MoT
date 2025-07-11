def solve_wuxing_factorial():
    """
    This function generates the C code for the Wuxing computer to calculate 100!
    and prints it in the required z[C] format.
    """
    memory_size = 170  # Total memory in 'D' for variables

    # C code for the Wuxing computer
    # The data types are chosen to minimize memory usage based on the problem description.
    # result[158]: 158 * 1D = 158D
    # i: 3D (char, for 2..100)
    # size: 3D (char, up to 158)
    # j: 3D (char, index up to 157)
    # carry: 3D (char, can be up to 999)
    # Total = 158 + 3 + 3 + 3 + 3 = 170D
    c_code = """#include <stdio.h>

void main() {
    // Defines variables using smallest possible memory footprint.
    digit result[158] = {0};
    char i;
    char size = 1;
    char j;
    char carry;

    result[0] = 1; // Initialize result to 1

    // Loop from i = 2 to 100, multiplying it with the result array.
    for (i = 2; i <= 100; i++) {
        carry = 0;
        // Multiply each digit of the result array with i.
        for (j = 0; j < size; j++) {
            carry += result[j] * i;
            result[j] = carry % 10;
            carry /= 10;
        }

        // Propagate the remaining carry, extending the array size if needed.
        while (carry > 0) {
            result[size] = carry % 10;
            carry /= 10;
            size++;
        }
    }

    // Print the final equation result: 100! = [the number]
    // The loop prints each digit, which is a number in the final result.
    printf("100! = ");
    
    // Print digits in reverse order. Use a while loop as 'j' is unsigned.
    j = size - 1;
    while (1) {
        printf("%d", result[j]);
        if (j == 0) {
            break;
        }
        j--;
    }
    printf("\\n");
}"""

    # Format and print the final answer as z[C]
    print(f"{memory_size}[{c_code}]")

solve_wuxing_factorial()