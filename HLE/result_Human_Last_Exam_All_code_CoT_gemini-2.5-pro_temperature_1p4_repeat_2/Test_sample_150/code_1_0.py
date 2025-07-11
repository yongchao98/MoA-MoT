import sys

def solve():
    """
    This function generates and prints an efficient C code solution for the Bit++ problem.
    """
    c_code = """#include <stdio.h>

/*
 * This is an efficient C program to solve the Bit++ problem.
 * It is designed to be fast and use minimal memory.
 *
 * Problem Analysis:
 * The initial Java program is incorrect because it only handles cases like "++X"
 * and fails on "X++" or "X--". This C code correctly handles all four cases.
 *
 * Efficiency Plan:
 * 1. Data Types: We use 'short' for the counter 'n' and result 'x' since their
 *    values are small (max 150), saving memory over 'int'.
 * 2. Logic: For any statement ("++X", "X++", "--X", "X--"), the character
 *    at index 1 is uniquely '+' or '-'. We use this for a simple check.
 * 3. Speed: To minimize character comparisons and avoid branching, we use an
 *    arithmetic trick. In ASCII, '+' is 43 and '-' is 45.
 *    - If the char is '+': x += (44 - 43) which means x is incremented by 1.
 *    - If the char is '-': x += (44 - 45) which means x is decremented by 1.
 *    This avoids an if/else structure.
 */
int main() {
    short n;
    // Use %hd to read into a short integer
    scanf("%hd", &n);

    short x = 0;
    // A buffer of size 4 is needed for a 3-char string + null terminator
    char statement[4];

    // Loop n times, decrementing n until it is zero
    while (n-- > 0) {
        scanf("%s", statement);
        // Use the arithmetic trick based on the middle character's ASCII value
        x += (44 - statement[1]);
    }

    // Print the final result. %hd is the format specifier for short.
    printf("%hd\\n", x);

    return 0;
}
"""
    print(c_code)

solve()