def solve_brick_tiling():
    """
    Calculates the solution for the brick tiling problem as specified.

    The solution is composed of two parts:
    x: Memory used (in bytes) by an optimal C program for the variables.
    o: The maximal covered area for the given N and M.
    """

    # The dimensions of the rectangle as per the problem statement.
    n = 1000000001
    m = 1000000001

    # Part 1: Calculate 'x', the memory usage.
    # An optimal C program for N=1000000001 would use 64-bit integers
    # to avoid overflow, as N*M is ~10^18.
    # A 64-bit integer (long long int in C) uses 8 bytes.
    # We need variables for n, m, and the output.
    bytes_per_variable = 8
    num_variables = 3  # for n, m, and the output
    memory_used_x = bytes_per_variable * num_variables

    # Part 2: Calculate 'o', the maximal covered area.
    # The area of a brick is 2. The maximum area must be a multiple of 2.
    # This can be calculated with integer division as (N * M) // 2 * 2.
    maximal_area_o = (n * m) // 2 * 2
    
    # Print the final result in the format "x:o"
    print(f"{memory_used_x}:{maximal_area_o}")

solve_brick_tiling()