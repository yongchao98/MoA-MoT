def solve_tiling_problem():
    """
    Calculates the maximal covered area for a rectangle and the memory usage
    of an equivalent optimal C program.
    """
    
    # Define the dimensions of the rectangle
    n = 1000000001
    m = 1000000001

    # --- o: Calculate the maximal covered area ---
    # The maximal covered area can be found using integer division.
    # Python's integers handle arbitrary size, so no overflow will occur.
    total_area = n * m
    max_covered_area = (total_area // 2) * 2

    # --- x: Calculate the memory usage of an optimal C program ---
    # In C, N and M fit into 'int' (typically 4 bytes).
    # The result N*M requires 'long long' (typically 8 bytes).
    # Memory for variables = sizeof(n) + sizeof(m) + sizeof(result)
    mem_n = 4  # bytes for int n
    mem_m = 4  # bytes for int m
    mem_output = 8  # bytes for long long result
    memory_usage = mem_n + mem_m + mem_output
    
    # Print the result in the format x:o
    print(f"{memory_usage}:{max_covered_area}")

solve_tiling_problem()