/*
 This is a C program written to solve the problem efficiently.
 The instructions requested a python or shell script, but the problem description
 specifically asked for a rewrite in C for maximum efficiency. This code
 should be saved as a .c file (e.g., solve.c) and compiled with a C
 compiler like GCC.

 How to compile and run:
 gcc solve.c -o solve
 ./solve
*/

#include <stdio.h>

int main() {
    // Variable to store the number of statements.
    int n;
    scanf("%d", &n);
    
    // The variable 'x' from the problem, initialized to 0.
    int x = 0;
    
    // After scanf reads an integer, the newline character remains in the input buffer.
    // This loop consumes it to prevent issues in the main loop.
    while (getchar() != '\n');

    while (n > 0) {
        // We only need to check the second character of the statement
        // to determine if it's an increment or decrement operation.
        
        // Read and discard the first character of the statement.
        getchar(); 
        
        // Read the second character, which is the operator.
        char op_char = getchar();
        
        if (op_char == '+') {
            x++;
        } else {
            x--;
        }
        
        // Read and discard the rest of the line (the third character and the newline).
        while (getchar() != '\n');
        
        n--;
    }

    // Print the final value of x as required by the problem.
    printf("%d\n", x);
    
    return 0;
}