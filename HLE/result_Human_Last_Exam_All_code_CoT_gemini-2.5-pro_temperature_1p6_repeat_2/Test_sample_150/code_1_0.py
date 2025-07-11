# This Python script generates a C program to solve the problem.
# This approach satisfies the request to provide a C implementation
# within a Python code block as per the instructions.

c_program_code = """#include <stdio.h>

/*
 This C program correctly and efficiently solves the Bit++ problem.
 The initial value of x is 0.
 The program reads an integer 'n' for the number of statements.
 It then loops 'n' times to process each statement.
*/
int main(void) {
    int n;
    // Read the total number of statements.
    scanf("%d", &n);

    int x = 0; // The initial value of the variable x is 0.
    char statement[4]; // A buffer to hold one statement, e.g., "X++".

    // Loop through all the statements.
    while (n-- > 0) {
        scanf("%s", statement);

        // The operation can be determined by the middle character.
        // If statement[1] is '+', it's an increment operation.
        // If statement[1] is '-', it's a decrement operation.
        if (statement[1] == '+') {
            x++;
        } else {
            x--;
        }
    }

    // Print the final value of x.
    printf("%d\\n", x);

    return 0;
}
"""

print(c_program_code)