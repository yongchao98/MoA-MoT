def solve_brick_problem():
    """
    Calculates the solution for the brick tiling problem.

    The function determines:
    1. The memory required for an optimal C program to solve the problem.
    2. The maximal area that can be covered on an NxM rectangle with 2x1 bricks.

    It then prints the result in the format 'memory:area'.
    """

    # The values for N and M as specified in the problem.
    n = 1000000001
    m = 1000000001

    # Step 1: Calculate the memory usage for variables in a C program.
    # The values of N, M, and their product (N*M) are large (~10^9 and ~10^18).
    # In C, these require a 64-bit integer type, 'long long', which uses 8 bytes.
    # An optimal program would have variables for n, m, and the result.
    # Memory for n (long long) = 8 bytes
    # Memory for m (long long) = 8 bytes
    # Memory for the output (long long) = 8 bytes
    memory_used_in_bytes = 8 + 8 + 8

    # Step 2: Calculate the maximal coverable area.
    # The area of one brick is 2x1 = 2.
    # Therefore, the total area covered must be an even number.
    # The total area of the rectangle is N * M.
    # The maximal coverable area is the largest even number less than or equal to N * M.
    # Using integer arithmetic, this is calculated as (N * M) // 2 * 2.
    # Python integers handle arbitrary size, so overflow is not an issue here.
    total_area = n * m
    maximal_covered_area = (total_area // 2) * 2

    # Step 3: Print the result in the required format "x:o".
    # The numbers in the final output are memory_used_in_bytes and maximal_covered_area.
    print(f"{memory_used_in_bytes}:{maximal_covered_area}")

solve_brick_problem()