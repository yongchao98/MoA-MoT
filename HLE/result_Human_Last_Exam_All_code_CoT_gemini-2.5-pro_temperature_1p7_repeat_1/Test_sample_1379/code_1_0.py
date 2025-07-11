def solve_brick_problem():
    """
    Calculates the solution for the brick tiling problem.

    The problem asks for two things:
    1. The memory used by an optimal C program for its variables (n, m, output).
    2. The maximal area that can be covered on an NxM grid with 2x1 bricks.

    N = M = 1000000001
    """

    # 1. Memory Calculation (for a hypothetical C program)
    # The values for N, M, and the result (N*M) are large and require
    # a 64-bit integer type in C, which is 'long long'.
    # A 'long long' variable uses 8 bytes of memory.
    # The C program would need variables for n, m, and the output.
    size_of_long_long_bytes = 8
    num_variables = 3  # n, m, output
    memory_x = size_of_long_long_bytes * num_variables

    # 2. Maximal Area Calculation
    n = 1000000001
    m = 1000000001

    # The area of one brick is 2. The total covered area must be even.
    # The maximum covered area is the largest even number less than or equal to the total area (n * m).
    # This can be calculated efficiently using integer division.
    total_area = n * m
    max_covered_area_o = (total_area // 2) * 2

    # The final equation to get the output `o` is (1000000001 * 1000000001) // 2 * 2.
    # Here we print each component of that equation before the final result.
    print(f"n = {n}")
    print(f"m = {m}")
    print(f"Final equation for area: ({n} * {m}) // 2 * 2")
    
    # Print the final result in the format x:o
    print(f"Final Answer: {memory_x}:{max_covered_area_o}")


solve_brick_problem()