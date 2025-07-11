def solve():
    """
    This script calculates the memory usage for an efficient C interpreter
    for the X++ language based on the provided constraints.
    """

    # Step 1 & 2: Determine the maximum number of statements (n).
    # The total number of characters in a program file must be <= 366.
    # The file format is:
    #   - The number of statements, n, as a string. Let its length be L.
    #   - A newline character (1 char).
    #   - n statements, each being 3 characters (e.g., "X++") plus a newline (1 char), totaling 4 chars per statement.
    #
    # The formula for total characters is: L + 1 + (n * 4) <= 366.
    #
    # Let's test for the number of digits in n (L):
    # - If n has 2 digits (L=2, e.g., n=10 to 99):
    #   2 + 1 + 4*n <= 366
    #   3 + 4*n <= 366
    #   4*n <= 363
    #   n <= 90.75
    #   So, the maximum 2-digit n is 90.
    #
    # - If n has 3 digits (L=3, e.g., n=100 to 999):
    #   3 + 1 + 4*n <= 366
    #   4 + 4*n <= 366
    #   4*n <= 362
    #   n <= 90.5
    #   This is a contradiction, as we assumed n is a 3-digit number (>=100).
    #
    # Therefore, the maximum number of statements is 90.
    max_n = 90

    # Step 3 & 4: Identify variables and select optimal data types for the C interpreter.
    # The efficient interpreter will read character by character, not line by line.

    # Variable 1: To store the number of statements, 'n'.
    # Range: 0 to 90.
    # Smallest type that can hold this is int8_t (1 byte, range -128 to 127).
    n_memory = 1

    # Variable 2: To store the result, 'x'.
    # Initial value is 0. It is incremented or decremented 'n' times.
    # Range: -90 to +90.
    # Smallest type that can hold this is also int8_t (1 byte, range -128 to 127).
    x_memory = 1

    # Variable 3: To store the character read from input, 'c'.
    # The interpreter reads one character at a time using getchar().
    # Since the problem specifies separate eoln() and eof() functions, we don't need
    # to check the return value of getchar() for EOF. A standard 'char' type is sufficient.
    # A 'char' in C uses 1 byte.
    c_memory = 1

    # Step 5: Calculate the total memory.
    total_memory = n_memory + x_memory + c_memory
    
    print("To write the most memory-efficient C interpreter, we need to determine the smallest possible data types for our variables.")
    print("\nStep 1: Find the maximum number of statements (n).")
    print(f"The program tape is limited to 366 characters. A program with 'n' statements and 'L' digits in 'n' has a length of L + 1 + (n * 4).")
    print(f"By solving '2 + 1 + 4*n <= 366', we find the maximum n is {max_n}.")

    print("\nStep 2: Choose data types for the variables.")
    print(f" - Variable 'n' (statement count): Max value is {max_n}. This fits in an int8_t, which is {n_memory} byte.")
    print(f" - Variable 'x' (the result): Value ranges from -{max_n} to +{max_n}. This also fits in an int8_t, which is {x_memory} byte.")
    print(f" - Variable 'c' (for getchar()): To read one character at a time, a char type is sufficient. This is {c_memory} byte.")
    
    print("\nStep 3: Calculate the total memory usage.")
    print("The total memory is the sum of the memory for each variable.")
    print(f"Memory for 'n' ({n_memory} byte) + Memory for 'x' ({x_memory} byte) + Memory for 'c' ({c_memory} byte) = {total_memory} bytes.")

solve()