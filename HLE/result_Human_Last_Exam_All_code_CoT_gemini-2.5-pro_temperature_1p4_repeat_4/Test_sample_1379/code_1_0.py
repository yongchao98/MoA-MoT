def solve_brick_problem():
    """
    This function calculates the memory usage and maximal covered area
    for the given problem and prints the result in the specified format.
    """
    
    # Input values for N and M
    n = 1000000001
    m = 1000000001

    # In an optimal C program, we use the smallest sufficient data types.
    # N and M fit into a 32-bit integer (4 bytes).
    bytes_for_n = 4
    bytes_for_m = 4

    # The product N*M is ~10^18, which requires a 64-bit integer (8 bytes).
    # The output variable must therefore be 64-bit.
    bytes_for_output = 8
    
    # The total memory usage for the specified variables.
    total_memory_bytes = bytes_for_n + bytes_for_m + bytes_for_output

    # The maximal covered area is the largest even number <= N*M.
    # This is calculated using integer division.
    total_area = n * m
    max_covered_area = (total_area // 2) * 2

    # Print the final result in the format x:o
    print(f"{total_memory_bytes}:{max_covered_area}")

solve_brick_problem()