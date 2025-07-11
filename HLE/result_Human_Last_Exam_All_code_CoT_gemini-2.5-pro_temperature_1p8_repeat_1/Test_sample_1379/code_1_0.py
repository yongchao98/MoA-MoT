def solve_brick_tiling():
    """
    Calculates the maximal area of an NxM rectangle that can be covered by 2x1 bricks,
    and the memory required for variables in an equivalent C program.

    The problem is solved for N = M = 1000000001.
    """
    # The dimensions of the rectangle.
    N = 1000000001
    M = 1000000001

    # In a C program for these values, N, M, and their product would require a
    # 64-bit integer type (long long int) to avoid overflow.
    # A long long int is 8 bytes. The memory for the three key variables is calculated.
    memory_for_n = 8  # bytes for a long long int
    memory_for_m = 8  # bytes for a long long int
    memory_for_output = 8 # bytes for a long long int
    total_memory_x = memory_for_n + memory_for_m + memory_for_output

    # The area of a single brick is 2x1=2.
    # Therefore, the maximal area that can be covered is the total area of the
    # rectangle rounded down to the nearest even number.
    # The equation for the output is: (N * M) // 2 * 2
    # In this final equation, you can see all the numbers involved: N, M, and 2.
    total_area = N * M
    maximal_covered_area_o = (total_area // 2) * 2

    # Print the final result in the format x:o
    print(f"{total_memory_x}:{maximal_covered_area_o}")

solve_brick_tiling()