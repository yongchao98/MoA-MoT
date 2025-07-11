def solve_tiling_problem():
    """
    Calculates the memory usage for an optimal C program and the maximal
    coverable area of a rectangle with 2x1 bricks.
    """
    # Given dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # Step 1: Calculate the memory usage for an optimal C program.
    # The values of N, M, and their product require a 64-bit integer type
    # to avoid overflow. In C, this is 'long long'.
    # A 'long long' typically uses 8 bytes of memory.
    # We need variables for n, m, and the output.
    size_of_long_long = 8  # in bytes
    memory_for_n = size_of_long_long
    memory_for_m = size_of_long_long
    memory_for_output = size_of_long_long
    total_memory_usage = memory_for_n + memory_for_m + memory_for_output

    # Step 2: Calculate the maximal coverable area.
    # The number of 2x1 bricks that can fit is floor(N * M / 2).
    # The maximal area is 2 * floor(N * M / 2).
    # In integer arithmetic (like Python's //), this is (n * m // 2) * 2.
    total_area = n * m
    # Using integer division // which is equivalent to floor division
    num_bricks = total_area // 2
    max_covered_area = num_bricks * 2

    # Step 3: Print the result in the specified format "memory:output".
    # The instruction "output each number in the final equation" is interpreted
    # as calculating the components of the final output string rather than
    # hardcoding them. The numbers in the final output are the total memory
    # usage and the maximum covered area.
    print(f"{total_memory_usage}:{max_covered_area}")

solve_tiling_problem()