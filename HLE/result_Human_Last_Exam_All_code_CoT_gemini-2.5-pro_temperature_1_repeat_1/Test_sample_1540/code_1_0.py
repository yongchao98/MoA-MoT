# The user's task is to write a memory-efficient C program.
# Per the instructions, I am providing the code inside a Python block.
# You can save the C code below into a file (e.g., compare.c),
# compile it with a C compiler (e.g., gcc compare.c -o compare),
# and then execute it (e.g., ./compare).

#include <stdio.h>
#include <ctype.h>

/*
 * This program compares two strings lexicographically, ignoring case.
 * It is designed to be memory-efficient by only storing the first
 * string in a buffer and reading the second string character-by-character.
 */
int main(void) {
    // A buffer for the first string. Max length 100 + 1 for the null terminator.
    char s1[101];

    // A loop counter. On most contest systems, an int is 4 bytes.
    int i = 0;

    // A variable to hold characters from the second string.
    // It is an 'int' to correctly handle the EOF value returned by getchar().
    int c2;

    // Read the first string from standard input.
    // scanf reads up to the first whitespace, leaving the newline in the buffer.
    scanf("%s", s1);

    // We must consume the newline character that separates the two input strings.
    getchar();

    // Loop until the end of the first string (indicated by the null terminator).
    while (s1[i]) {
        // Read the next character from the input stream (which is the second string).
        c2 = getchar();

        // Compare the lowercase versions of the characters.
        if (tolower(s1[i]) < tolower(c2)) {
            printf("-1\n");
            return 0; // Exit after finding the result.
        }
        if (tolower(s1[i]) > tolower(c2)) {
            printf("1\n");
            return 0; // Exit after finding the result.
        }

        // If characters are equal, move to the next one.
        i++;
    }

    // If the loop finishes, all corresponding characters were equal.
    // The problem guarantees the strings are of the same length.
    printf("0\n");

    return 0;
}
<<<109>>>