def solve_brick_problem():
    """
    Calculates the memory usage for a C program and the maximal tileable area.
    """
    
    # Step 1: Define the dimensions of the rectangle as per the problem.
    # N = 1,000,000,001
    # M = 1,000,000,001
    n = 1000000001
    m = 1000000001

    # Step 2: Calculate the memory required for an optimal C program.
    # To store numbers of this magnitude (N*M is ~10^18), an `unsigned long long`
    # is required in C, which is 8 bytes.
    # The memory is for three variables: n, m, and the output.
    memory_for_n = 8
    memory_for_m = 8
    memory_for_output = 8
    total_memory_bytes = memory_for_n + memory_for_m + memory_for_output
    
    # Step 3: Calculate the maximal covered area.
    # The area of one brick is 2. The total area that can be covered must
    # be an even number. The maximal area is the largest even number
    # less than or equal to the total area N * M.
    # Using integer division, this is calculated as (N * M) // 2 * 2.
    # Python's integers handle arbitrary precision, so overflow is not an issue here.
    total_area = n * m
    maximal_covered_area = (total_area // 2) * 2

    # Step 4: Print the result in the specified "x:o" format.
    print(f"{total_memory_bytes}:{maximal_covered_area}")

solve_brick_problem()