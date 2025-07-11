def solve():
    """
    This function generates the final answer string in the format z[C].
    z: The smallest memory size in decimal digits (D).
    C: The C program code to calculate 100!.
    """
    memory_size = "172D"

    c_code = """#include <stdio.h>

void main() {
    // Total variable memory: 158D (result) + 3D (i) + 3D (len) 
    // + 3D (j) + 3D (product) + 2D (carry) = 172D
    digit result[158] = {0};
    char i, len, j, product;
    short carry;

    // Initialize result to 1
    result[0] = 1;
    len = 1;

    // Multiply by i from 2 to 100
    for (i = 2; i <= 100; i++) {
        carry = 0;
        for (j = 0; j < len; j++) {
            product = result[j] * i + carry;
            result[j] = product % 10;
            carry = product / 10;
        }
        
        // Handle remaining carry by extending the number
        while (carry > 0) {
            result[len] = carry % 10;
            carry /= 10;
            len++;
        }
    }

    // Print each digit of the final number, which is stored in reverse
    for (j = len - 1; j >= 0; j--) {
        printf("%d", result[j]);
    }
}"""

    # The prompt asks for the answer in the format z[C]
    final_answer = f"{memory_size}[{c_code}]"
    print(final_answer)

solve()