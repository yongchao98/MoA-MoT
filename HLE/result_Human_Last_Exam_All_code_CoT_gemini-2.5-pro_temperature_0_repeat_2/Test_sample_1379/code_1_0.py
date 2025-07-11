def solve_brick_problem():
    """
    Calculates the memory usage for an optimal C program and the maximal
    coverable area for a given rectangle, then prints the result.
    """
    # Given dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # --- Step 1: Calculate memory usage (x) ---
    # An optimal C program must use a 64-bit integer type ('long long') to
    # store the product of n and m without overflow. A 'long long' is 8 bytes.
    # We account for memory for n, m, and the output variable.
    size_of_long_long_in_bytes = 8
    memory_for_n = size_of_long_long_in_bytes
    memory_for_m = size_of_long_long_in_bytes
    memory_for_output = size_of_long_long_in_bytes
    x = memory_for_n + memory_for_m + memory_for_output

    # --- Step 2: Calculate maximal area (o) ---
    # The maximal area is floor((n * m) / 2) * 2.
    # Using integer division, this is (n * m) // 2 * 2.
    total_area = n * m
    o = (total_area // 2) * 2

    # --- Step 3: Output the result in the format x:o ---
    # The instruction "output each number in the final equation" is interpreted
    # as printing the calculated values of x and o in the final string.
    print(f"{x}:{o}")

solve_brick_problem()