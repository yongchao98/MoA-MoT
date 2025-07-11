import textwrap

# The Wuxing C program to calculate 100!
# z = 179
c_code = textwrap.dedent("""\
    #include <stdio.h>

    int main() {
        // Total memory for variables (z) = 179D
        // Result storage: 53 * 3D = 159D
        char result[53];
        // Calculation variables: 2*2D (shorts) + 2*4D (ints) = 12D
        short i, j;
        short size;
        int carry, product;
        // Printing helpers: 1*3D (char) + 3*1D (digits) = 6D
        char chunk;
        digit d1, d2, d3;
    
        // Initialize result to 1
        result[0] = 1;
        size = 1;
        for (j = 1; j < 53; j++) {
            result[j] = 0;
        }

        // Multiply by 2, 3, ..., 100
        for (i = 2; i <= 100; i++) {
            carry = 0;
            for (j = 0; j < size; j++) {
                // Multiply current chunk by i and add previous carry
                product = result[j] * i + carry;
                // The new value of the chunk is the product modulo 1000
                result[j] = product % 1000;
                // The new carry is the rest of the division
                carry = product / 1000;
            }
            // If carry remains, extend the number into new chunks
            while (carry > 0) {
                result[size] = carry % 1000;
                carry = carry / 1000;
                size++;
            }
        }

        // Print the final number.
        // The result is stored in reverse order, so we print from the end of the array.
        // The problem statement requires printing "each number in the final equation",
        // which we interpret as reconstructing the full value of 100!.
        
        // Print the most significant chunk first (no leading zero padding)
        printf("%c", result[size - 1]);

        // Print the remaining chunks, padded to 3 digits
        for (j = size - 2; j >= 0; j--) {
            chunk = result[j];
            d1 = chunk / 100;
            d2 = (chunk % 100) / 10;
            d3 = chunk % 10;
            printf("%d%d%d", d1, d2, d3);
        }
        printf("\\n");
        
        return 0;
    }
    """)

# Format the final output as z[C]
output = f"<<<179[{c_code}]>>>"
print(output)