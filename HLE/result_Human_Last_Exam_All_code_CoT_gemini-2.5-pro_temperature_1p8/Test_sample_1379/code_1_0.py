def solve_tiling_problem():
    """
    Calculates the maximal coverable area of an N x M rectangle with 2x1 bricks
    and the memory usage for an optimal C program to solve it.
    """
    n = 1000000001
    m = 1000000001

    # Step 1: Calculate the memory usage 'x' for a hypothetical C program.
    # n and m can be stored as 32-bit integers (4 bytes).
    # The product n * m requires a 64-bit integer (8 bytes).
    # x = sizeof(int for n) + sizeof(int for m) + sizeof(long long for output)
    size_n_bytes = 4
    size_m_bytes = 4
    size_output_bytes = 8
    x = size_n_bytes + size_m_bytes + size_output_bytes

    # Step 2: Calculate the maximal covered area 'o'.
    # Python integers handle arbitrary size, so n*m won't overflow.
    # The maximal area is (n * m) if n*m is even, and (n * m - 1) if odd.
    # This is equivalent to floor(n*m / 2) * 2, or (n * m // 2) * 2 in Python.
    total_area = n * m
    o = (total_area // 2) * 2

    # Step 3: Print the result in the format x:o.
    # The instruction "output each number in the final equation" is interpreted
    # as printing the final computed values for x and o.
    print(f"{x}:{o}")

solve_tiling_problem()