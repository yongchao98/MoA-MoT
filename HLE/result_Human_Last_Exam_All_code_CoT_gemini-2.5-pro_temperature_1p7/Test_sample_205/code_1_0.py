import sys

def solve():
    """
    This function generates the C code for the Wuxing computer to calculate 100!
    and prints it along with its memory usage in the required z[C] format.
    """

    # z: Calculated memory usage in decimal digits (D)
    # result[160] -> 160*1D = 160D
    # i, j, size, carry (all char) -> 4 * 3D = 12D
    # Total z = 160 + 12 = 172D
    z = 172

    # C: The C code for the Wuxing architecture.
    # It is formatted for readability.
    C_code = r'''void main() {
    // Variable declarations. Types are chosen for minimal memory footprint.
    // result: stores factorial digits. 160 is safe for 100! (~158 digits). Memory: 160 * 1D = 160D
    // i, j, size, carry: counters and temporary variables. 
    // `char` (0-999) is sufficient for all of them. Memory: 4 * 3D = 12D
    // Total memory `z` = 160 + 12 = 172 D.
    digit result[160];
    char i, j, size, carry;

    // Initialize result to 1.
    result[0] = 1;
    size = 1;

    // Set other digits to 0. We can reuse 'j' as the counter.
    for (j = 1; j < 160; j = j + 1) {
        result[j] = 0;
    }

    // Calculate 100! by multiplying with numbers from 2 to 100.
    for (i = 2; i <= 100; i = i + 1) {
        carry = 0;
        // Multiply each digit of the current result by 'i'.
        for (j = 0; j < size; j = j + 1) {
            // Product of a digit * i + previous carry will fit in 'char' (0-999)
            // e.g., 9 * 100 + carry < 999
            carry = result[j] * i + carry;
            result[j] = carry % 10;
            carry = carry / 10;
        }

        // If there is a remaining carry, extend the number's size.
        while (carry > 0) {
            result[size] = carry % 10;
            size = size + 1;
            carry = carry / 10;
        }
    }
    
    // Print the final result in the format "100! = [number]".
    printf("100! = ");
    
    // Print digits from most to least significant.
    // A 'while' loop is safer for a potentially unsigned counter 'j'.
    j = size;
    while (j > 0) {
        j = j - 1;
        printf("%d", result[j]);
    }
    printf("\n");
}'''

    # Print the final answer in the specified z[C] format.
    # We use a file-like object to prevent any extra newlines if the
    # code string has one at the end.
    sys.stdout.write(f"{z}[{C_code}]")

solve()