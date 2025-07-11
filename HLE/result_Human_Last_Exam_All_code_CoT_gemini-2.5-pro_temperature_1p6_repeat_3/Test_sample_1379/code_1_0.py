def solve_brick_problem():
    """
    Calculates the memory usage for an optimal C program and the maximal
    coverable area for a given rectangle, then prints the result.
    """
    # The dimensions of the rectangle
    N = 1000000001
    M = 1000000001

    # --- Memory Calculation (x) ---
    # An optimal C program for these inputs needs 64-bit integers (`long long`)
    # to store the product N * M, which is ~10^18.
    # sizeof(long long) is typically 8 bytes on modern systems.
    # The variables are n, m, and the output.
    bytes_for_n = 8
    bytes_for_m = 8
    bytes_for_output = 8
    memory_x = bytes_for_n + bytes_for_m + bytes_for_output

    # --- Maximal Area Calculation (o) ---
    # The maximal area that can be covered by 2x1 bricks is the largest
    # even number less than or equal to the total area (N * M).
    # This can be calculated using integer division.
    # The equation to calculate the output is:
    total_area = N * M
    output_o = (total_area // 2) * 2

    # Print the final result in the format x:o
    print(f"{memory_x}:{output_o}")

solve_brick_problem()