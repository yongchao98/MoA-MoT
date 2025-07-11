def solve_brick_problem():
    """
    Calculates the maximal area coverable by 2x1 bricks on an NxM rectangle
    and the memory usage of an optimal C program to solve it.
    """

    # Given values for N and M
    N = 1000000001
    M = 1000000001

    # --- Part 1: Calculate the maximal covered area (o) ---
    # The total area is N * M.
    # The maximal covered area is the largest even number <= N * M.
    # Python's integers handle arbitrary size, so we can calculate directly.
    # Using integer division `//` achieves this: (N * M // 2) * 2
    output_o = (N * M) // 2 * 2

    # --- Part 2: Calculate memory usage for C variables (x) ---
    # An optimal C program must handle the given values and the result.
    # n = 1,000,000,001. A 'long long' (8 bytes) is a robust choice for inputs.
    # m = 1,000,000,001. A 'long long' (8 bytes) is also used for m.
    # output = ~10^18. This requires a 'long long' (8 bytes).
    # Total memory = sizeof(n) + sizeof(m) + sizeof(output)
    size_of_long_long = 8
    memory_x = size_of_long_long + size_of_long_long + size_of_long_long
    
    # --- Final Output ---
    # Print the result in the specified "x:o" format.
    # The numbers in the final output are the calculated memory_x and output_o.
    print(f"{memory_x}:{output_o}")

solve_brick_problem()