import math

def generate_wuxing_c_code():
    """
    This function generates the final output string in the format z[C].
    z = total memory used by variables in Decimal Digits (D).
    C = the C code for the Wuxing architecture.
    """

    # Calculate 100! to include it as a comment in the C code.
    factorial_100 = math.factorial(100)
    
    # Calculate the total memory size 'z'.
    # result[158] -> 158 * 1D = 158D
    # 5 char variables (result_size, i, j, carry, product) -> 5 * 3D = 15D
    # Total z = 158 + 15 = 173D
    z = 173

    # The C program code 'C'.
    # Note: Backslashes in printf must be escaped (e.g., "\\n") for the string.
    c_code = f"""#include <stdio.h>

void main() {{
    // This program calculates 100! using an array to handle the large number.
    // The final result is:
    // {factorial_100}

    digit result[158]; // Array to store the 158 digits of the result.
    char result_size;  // Tracks the current number of digits in the result.
    char i, j;         // Loop counters.
    char carry;        // Stores the carry-over during multiplication.
    char product;      // Stores the intermediate product of a digit multiplication.

    // Initialize the result to 1.
    result_size = 1;
    result[0] = 1;

    // Main loop: multiply the current result by i, from 2 to 100.
    for (i = 2; i <= 100; i++) {{
        carry = 0;
        for (j = 0; j < result_size; j++) {{
            // Multiply the current digit by i and add the previous carry.
            product = result[j] * i + carry;
            // The new digit is the remainder of the product divided by 10.
            result[j] = product % 10;
            // The new carry is the integer division of the product by 10.
            carry = product / 10;
        }}

        // If there is a carry left over, append its digits to the result array.
        while (carry > 0) {{
            result[result_size] = carry % 10;
            carry = carry / 10;
            result_size++;
        }}
    }}

    // Print the final equation as "100! = [number]".
    printf("100! = ");
    
    // Print the digits from the result array in reverse order.
    // The loop counter 'j' is reused for this purpose.
    j = result_size;
    while (j > 0) {{
        j--;
        printf("%d", result[j]);
    }}
    printf("\\n");
}}"""

    # Format the final answer string as z[C].
    final_answer = f"{z}[{c_code}]"
    print(final_answer)

generate_wuxing_c_code()