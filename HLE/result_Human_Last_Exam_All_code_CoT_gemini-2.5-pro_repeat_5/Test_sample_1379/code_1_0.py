def solve_tiling_problem():
    """
    Calculates the maximal area of an NxM rectangle that can be covered by 2x1 bricks
    and the memory required for an optimal C program to store the variables.
    """
    
    # Given values for N and M
    n = 1000000001
    m = 1000000001

    # Step 1: Calculate the maximal covered area.
    # The number of bricks is floor((n * m) / 2).
    # The area covered is 2 * floor((n * m) / 2).
    # In Python, this is equivalent to integer division.
    total_area = n * m
    max_covered_area = (total_area // 2) * 2

    # Step 2: Calculate memory usage for an equivalent optimal C program.
    # n = 1,000,000,001 fits in a 4-byte unsigned int.
    mem_n_bytes = 4
    # m = 1,000,000,001 fits in a 4-byte unsigned int.
    mem_m_bytes = 4
    # The output, ~10^18, requires an 8-byte unsigned long long.
    mem_output_bytes = 8
    
    total_memory_bytes = mem_n_bytes + mem_m_bytes + mem_output_bytes
    
    # Print the result in the format x:o
    print(f"{total_memory_bytes}:{max_covered_area}")

solve_tiling_problem()