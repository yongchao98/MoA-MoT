def solve():
    """
    This function generates the C code for calculating 100! on the WXVM
    and prints it along with the calculated memory usage.
    """

    # Memory usage calculation based on Wuxing C data types:
    # - The result of 100! has 158 digits. An array of size 160 is safe.
    #   digit result[160]; -> 160 * 1D = 160D
    # - Loop counter 'i' (2 to 100) fits in a 'char'.
    #   char i; -> 3D
    # - 'num_digits' tracks the number of digits (max ~158), fits in a 'char'.
    #   char num_digits; -> 3D
    # - Inner loop counter 'j' (max ~158) fits in a 'char'.
    #   char j; -> 3D
    # - 'carry' variable will temporarily hold the product (e.g., 9*100 + carry),
    #   which is at most 999. This fits in a 'char'.
    #   char carry; -> 3D
    # Total memory = 160 + 3 + 3 + 3 + 3 = 172D
    memory_size = 172

    # The C code for the Wuxing virtual machine.
    # Note: 'digit' is a 1D type (0-9), 'char' is a 3D unsigned type (0-999).
    c_code = """
void main() {
    // This program calculates 100! using an array of digits
    // to handle the large result. Memory usage is minimized by
    // selecting the smallest possible data types for all variables.

    // Variable Declarations (Total Memory: 172D)
    digit result[160]; // 160 * 1D = 160D
    char i;            // 3D (for loop 2-100)
    char num_digits;   // 3D (max ~158)
    char j;            // 3D (for inner loops)
    char carry;        // 3D (max product ~999)

    // Initialize result to 1
    result[0] = 1;
    num_digits = 1;

    // Calculate 100! by multiplying from 2 to 100
    for (i = 2; i <= 100; i++) {
        carry = 0;
        // Multiply the number in 'result' array by 'i'
        for (j = 0; j < num_digits; j++) {
            // 'carry' temporarily stores the product before being updated
            carry = result[j] * i + carry;
            result[j] = carry % 10;
            carry = carry / 10;
        }

        // Propagate the remaining carry to extend the number
        while (carry > 0) {
            result[num_digits] = carry % 10;
            num_digits++;
            carry = carry / 10;
        }
    }

    // Print the final equation: 100! = <result>
    printf("100! = ");
    
    // Print each digit of the result.
    // A while loop is used because 'char' is unsigned, and a
    // for loop like 'for(j=num_digits-1; j>=0; j--)' would be infinite.
    j = num_digits;
    while (j > 0) {
        j--;
        printf("%d", result[j]);
    }
    printf("\\n");
}"""

    # The final output format is z[C]
    # The strip() method is used to remove leading/trailing whitespace from the C code string.
    final_answer = f"<<<{memory_size}[{c_code.strip()}
]>>>"
    print(final_answer)

solve()