import sys

def calculate_memory_usage():
    """
    This function calculates the minimum memory required for the C interpreter's variables.
    """
    
    # Step 1: Determine the longest program.
    # A program is: n_string + newline + n * (statement + newline)
    # Total chars <= 366.
    # Statement length is 3. Newline is 1 char. Line length is 4.
    # Total chars = len(str(n)) + 1 + n * 4 <= 366
    # By testing values:
    # If n = 90, Total chars = len("90") + 1 + 90 * 4 = 2 + 1 + 360 = 363 <= 366. (OK)
    # If n = 91, Total chars = len("91") + 1 + 91 * 4 = 2 + 1 + 364 = 367 > 366. (Too long)
    # So, the maximum number of statements is n = 90.

    # Step 2: Identify variables for the most memory-efficient C interpreter.
    # The C interpreter needs variables to store:
    # - The number of statements, n.
    # - The value of the X++ variable, x.
    # - The character being read from the tape, c.
    
    # Step 3: Determine the minimum size for each variable in bytes.
    # Variable 'n' needs to store a maximum value of 90.
    # A signed 8-bit integer (int8, 1 byte) can store -128 to 127, which is sufficient.
    size_of_n_in_bytes = 1

    # Variable 'x' accumulates the result. Its value will be between -90 and 90.
    # A signed 8-bit integer (int8, 1 byte) is also sufficient for this range.
    size_of_x_in_bytes = 1

    # Variable 'c' is used to store the return value of getchar().
    # The getchar() function in C returns an 'int' to be able to represent all
    # characters plus the special End-Of-File (EOF) value.
    # The problem specifies that an 'int' in Xland is 32 bits.
    # 32 bits = 4 bytes.
    size_of_c_in_bytes = 4

    # Step 4: Calculate the total memory.
    total_memory = size_of_n_in_bytes + size_of_x_in_bytes + size_of_c_in_bytes

    # Print the detailed breakdown of the calculation.
    print(f"The minimum memory in bytes is estimated based on the C variables required:")
    print(f"1. Memory for variable 'n' (to store the number of statements, max 90): {size_of_n_in_bytes} byte.")
    print(f"2. Memory for variable 'x' (to store the result, range [-90, 90]): {size_of_x_in_bytes} byte.")
    print(f"3. Memory for variable 'c' (to read characters via getchar()): {size_of_c_in_bytes} bytes.")
    print("\nEquation for total memory usage:")
    print(f"{size_of_n_in_bytes} + {size_of_x_in_bytes} + {size_of_c_in_bytes} = {total_memory}")
    
    # Storing the final result in the requested format for the platform.
    # The problem asks for the answer to be returned at the end.
    sys.stdout.write(f"\n<<<{total_memory}>>>")


calculate_memory_usage()