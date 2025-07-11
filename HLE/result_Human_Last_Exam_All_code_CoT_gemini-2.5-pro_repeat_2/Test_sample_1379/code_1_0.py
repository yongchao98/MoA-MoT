def solve_tiling_problem():
    """
    Calculates the maximal coverable area of a rectangle with 2x1 bricks
    and the memory usage of an optimal C program to solve it.
    """
    
    # The dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # --- Calculate the maximal area (o) ---
    # Since both n and m are odd, the total area n*m is odd.
    # The maximal coverable area is n*m - 1.
    # Python's integers handle arbitrary size, so overflow is not an issue here.
    output_o = n * m - 1

    # --- Calculate the memory usage (x) ---
    # In an optimal C program, we would use the smallest possible data types.
    # n (1,000,000,001) fits in a 4-byte `int`.
    # m (1,000,000,001) fits in a 4-byte `int`.
    # The result `o` (~10^18) requires an 8-byte `long long`.
    # Total memory for variables = 4 (for n) + 4 (for m) + 8 (for the output).
    output_x = 4 + 4 + 8

    # Print the final result in the format x:o
    print(f"{output_x}:{output_o}")

solve_tiling_problem()