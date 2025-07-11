def solve_tiling_problem():
    """
    Calculates the maximal area that can be covered by 2x1 bricks on a NxM
    rectangle and estimates the memory usage of an equivalent C program.
    """
    # The dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # The area of a brick is 2. The total area of the rectangle is n * m.
    # The maximum number of bricks is floor((n * m) / 2).
    # The maximal covered area is floor((n * m) / 2) * 2.
    # In integer arithmetic, this is (n * m) // 2 * 2.
    total_area = n * m
    covered_area = (total_area // 2) * 2
    
    # Estimate memory usage for a C program.
    # To store n, m, and the result, we need to consider the data types.
    # n and m are ~10^9, which fit in a 32-bit signed int (4 bytes).
    # However, the product n * m is ~10^18, which requires a 64-bit integer
    # (long long in C, 8 bytes). To avoid overflow during calculation,
    # n and m would be cast to long long.
    # So, we'd have three 'long long' variables in a straightforward C program.
    # sizeof(long long) is typically 8 bytes.
    # Memory for n (long long) = 8 bytes
    # Memory for m (long long) = 8 bytes
    # Memory for output (long long) = 8 bytes
    # Total memory x:
    memory_x = 8 + 8 + 8
    
    # The output o:
    output_o = covered_area
    
    # Print the result in the format x:o
    print(f"{memory_x}:{output_o}")

solve_tiling_problem()