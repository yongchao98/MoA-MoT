def solve_brick_tiling():
    """
    Calculates the maximal coverable area of a rectangle with 2x1 bricks
    and the memory required for a C implementation.
    """
    # Given dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # Step 1: Calculate the maximal covered area 'o'.
    # The total area is n * m. Since each brick has an area of 2,
    # the maximum covered area is the largest even number <= n * m.
    # This can be calculated using integer division.
    total_area = n * m
    maximal_covered_area = (total_area // 2) * 2

    # Step 2: Calculate the memory usage 'x' for a hypothetical C program.
    # The values for n, m, and the result (~10^18) require a 64-bit integer.
    # In C, this corresponds to 'long long', which is typically 8 bytes.
    # Memory for n: 8 bytes
    # Memory for m: 8 bytes
    # Memory for the output variable: 8 bytes
    memory_for_n = 8
    memory_for_m = 8
    memory_for_output = 8
    total_memory_usage = memory_for_n + memory_for_m + memory_for_output
    
    # Step 3: Print the result in the format x:o.
    # The prompt requires outputting each number in the "final equation",
    # which is interpreted as the 'x:o' format.
    print(f"{total_memory_usage}:{maximal_covered_area}")

solve_brick_tiling()