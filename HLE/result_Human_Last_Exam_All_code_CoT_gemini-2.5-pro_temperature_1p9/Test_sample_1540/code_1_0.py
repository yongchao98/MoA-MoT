import sys

def generate_c_code_and_analysis():
    """
    This function generates and prints the most memory-efficient C program 
    for the described problem, along with its variable memory analysis.
    """

    # The memory-efficient C code, formatted as a string.
    c_code = r"""#include <stdio.h>
#include <ctype.h>

// This is the most memory-efficient C program to solve the problem.
int main() {
    // Each buffer stores up to 100 characters plus a null terminator.
    // sizeof(char) is 1 byte.
    // Memory usage for s1: 101 * 1 = 101 bytes.
    char s1[101];

    // Memory usage for s2: 101 * 1 = 101 bytes.
    char s2[101];

    // Read the two strings from the input.
    scanf("%s", s1);
    scanf("%s", s2);

    // For the loop index, an 'unsigned char' is used instead of an 'int'.
    // An unsigned char uses 1 byte and can hold values from 0-255,
    // which is sufficient since the string length is at most 100.
    // Memory usage for i: 1 byte.
    unsigned char i = 0;

    // Loop until the null terminator of the first string is reached.
    while (s1[i]) {
        // Use tolower() for case-insensitive comparison directly.
        // This avoids creating extra temporary variables.
        if (tolower(s1[i]) < tolower(s2[i])) {
            printf("-1\n");
            return 0;
        }
        if (tolower(s1[i]) > tolower(s2[i])) {
            printf("1\n");
            return 0;
        }
        i++;
    }

    // If the loop completes, the strings are equal.
    printf("0\n");

    return 0;
}
"""

    # --- Calculation details ---
    s1_memory = 101
    s2_memory = 101
    index_memory = 1
    total_memory = s1_memory + s2_memory + index_memory

    # --- Print the analysis to the user ---
    print("--- Requested C Program ---")
    print("Below is the most memory-efficient C code for the given task:")
    print("-------------------------------------------------------------")
    print(c_code)
    print("-------------------------------------------------------------\n")
    print("--- Variable Memory Usage Analysis (m) ---")
    print("The total memory 'm' in bytes used by the program's variables is:")
    print(f"* Variable `s1` (char array of size 101): {s1_memory} bytes")
    print(f"* Variable `s2` (char array of size 101): {s2_memory} bytes")
    print(f"* Variable `i` (loop index of type unsigned char): {index_memory} byte")
    print("\nThe final equation for total variable memory is:")
    print(f"m = {s1_memory} + {s2_memory} + {index_memory} = {total_memory} bytes")


if __name__ == '__main__':
    generate_c_code_and_analysis()
