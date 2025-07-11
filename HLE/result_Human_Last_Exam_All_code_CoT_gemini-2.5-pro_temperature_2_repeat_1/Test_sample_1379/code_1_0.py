def solve_brick_problem():
    """
    Calculates the maximal area of an NxM rectangle that can be covered by 2x1 bricks
    and the memory usage of an optimal C program to solve it.
    """
    # The dimensions of the rectangle as per the problem
    n = 1000000001
    m = 1000000001

    # --- Area Calculation ---
    # The total area of the rectangle
    total_area = n * m

    # Each brick has an area of 2. The maximum number of bricks is floor(total_area / 2).
    # The maximal covered area is this number of bricks multiplied by 2.
    # In integer arithmetic, this is (total_area // 2) * 2.
    max_covered_area = (total_area // 2) * 2

    # --- Memory Usage Analysis for a C Program ---
    # N = 1,000,000,001 fits in a 32-bit signed integer (int).
    # In most common systems (like x86-64 Linux/Windows), sizeof(int) is 4 bytes.
    memory_for_n_bytes = 4

    # M = 1,000,000,001 also fits in a 32-bit signed integer (int).
    # sizeof(int) is 4 bytes.
    memory_for_m_bytes = 4

    # The output (n * m) can be up to ~10^18. This requires a 64-bit integer type.
    # In C, this is `long long`. sizeof(long long) is 8 bytes.
    memory_for_output_bytes = 8

    # Total memory for the variables in an optimal C program.
    total_memory_bytes = memory_for_n_bytes + memory_for_m_bytes + memory_for_output_bytes
    
    # Print the final result in the format x:o
    # x = total_memory_bytes
    # o = max_covered_area
    print(f"{total_memory_bytes}:{max_covered_area}")

solve_brick_problem()