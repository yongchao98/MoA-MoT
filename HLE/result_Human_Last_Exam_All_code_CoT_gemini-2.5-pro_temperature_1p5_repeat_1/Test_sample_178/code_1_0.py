def solve():
    """
    This script calculates the memory usage for the most efficient C interpreter
    for the X++ language, based on the problem's constraints.
    """

    # Step 1: Determine the maximum number of statements (n).
    # The total number of characters read from the tape must be <= 366.
    # A program has n statements. Let d be the number of digits in n.
    # Total chars = (d for the number n + 1 for newline) + (n statements * 4 chars/statement)
    # d + 1 + 4n <= 366
    # Let's test values of n:
    # If n=90, d=2. Total = 2 + 1 + 4*90 = 3 + 360 = 363. (Valid)
    # If n=91, d=2. Total = 2 + 1 + 4*91 = 3 + 364 = 367. (Invalid)
    # So, the maximum number of statements is 90.
    max_n = 90

    # Step 2: Identify variables for the most efficient C interpreter and their sizes.
    # A memory-efficient C program needs the following variables:
    
    # - A variable 'x' to store the result.
    #   With max_n = 90, the value of x will be in the range [-90, 90].
    #   The smallest C type that can hold this is a signed 8-bit integer (int8_t).
    size_x = 1 # bytes for int8_t

    # - A variable 'n' to store the number of statements.
    #   The maximum value is 90.
    #   The smallest C type that can hold this is an unsigned 8-bit integer (uint8_t).
    #   This variable can also be used as the loop counter (e.g., while(n--)),
    #   avoiding the need for a separate loop counter variable.
    size_n = 1 # bytes for uint8_t

    # - A variable 'c' to read characters using getchar().
    #   The getchar() function in C returns an 'int' to handle all character
    #   values plus the special End-Of-File (EOF) marker.
    #   The problem states that an 'int' is 32 bits.
    size_c = 4 # bytes for a 32-bit int

    # Step 3: Calculate the total memory.
    total_memory = size_x + size_n + size_c

    # Step 4: Print the detailed explanation and calculation.
    print("Estimating memory usage for the most efficient C interpreter:")
    print(f"1. The maximum number of statements (n) possible in a program up to 366 characters is {max_n}.")
    print("\n2. The interpreter requires three variables with minimum possible sizes:")
    print(f"   - Variable 'x' (result): Range is [{-max_n}, {max_n}]. Needs a signed 8-bit integer (int8_t), which is {size_x} byte.")
    print(f"   - Variable 'n' (statement count): Max value is {max_n}. Needs an unsigned 8-bit integer (uint8_t), which is {size_n} byte.")
    print(f"   - Variable 'c' (for getchar()): Must be a 32-bit 'int' as specified to handle EOF, which is {size_c} bytes.")
    
    print("\n3. The total memory is the sum of the sizes of these variables.")
    print(f"Total memory = size('x') + size('n') + size('c')")
    print(f"Total memory = {size_x} + {size_n} + {size_c} = {total_memory} bytes.")

solve()