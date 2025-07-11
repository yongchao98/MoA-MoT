def solve_brick_problem():
    """
    Calculates the memory usage for an optimal C program and the maximal
    covered area for a given N x M rectangle.
    """
    # Define the dimensions of the rectangle
    N = 1000000001
    M = 1000000001

    # --- Memory Calculation (x) ---
    # An optimal C program would use 64-bit integers (long long) to avoid
    # overflow when calculating the area N*M. A long long typically uses 8 bytes.
    # We need variables for n, m, and the output area.
    bytes_per_variable = 8
    num_variables = 3  # n, m, area
    memory_used_bytes = bytes_per_variable * num_variables

    # --- Area Calculation (o) ---
    # The total area is N * M. Each brick covers 2 units.
    # The number of bricks is the integer division of the total area by 2.
    # The maximum covered area is the number of bricks * 2.
    # In Python, this is `(N * M) // 2 * 2`.
    # This is equivalent to N*M if N*M is even, and N*M - 1 if N*M is odd.
    # Since N and M are both odd, their product is odd.
    max_covered_area = (N * M // 2) * 2

    # --- Output ---
    # Print the result in the format x:o
    print(f"{memory_used_bytes}:{max_covered_area}")

solve_brick_problem()