def solve_tiling():
    """
    Calculates the maximal area of an NxM rectangle that can be covered by 2x1 bricks,
    and the memory required to store the variables in a corresponding C program.
    """
    # Given dimensions of the rectangle
    N = 1000000001
    M = 1000000001

    # Calculate the maximal covered area.
    # The area of each brick is 2. The total area that can be covered must be even.
    # The total area of the rectangle is N * M.
    # The number of bricks is floor((N * M) / 2).
    # So, the maximal covered area is floor((N * M) / 2) * 2.
    # Python's '//' operator performs integer division (floor division).
    total_area = N * M
    max_covered_area = (total_area // 2) * 2

    # Calculate the memory usage (x) for a corresponding C program.
    # N and M are > 2^31 - 1, so they require a 64-bit integer type (long long in C).
    # On most modern systems, sizeof(long long) is 8 bytes.
    memory_for_n_bytes = 8
    
    memory_for_m_bytes = 8

    # The result, N*M, is approx 10^18, which fits in an unsigned 64-bit integer.
    # sizeof(unsigned long long) is 8 bytes.
    memory_for_output_bytes = 8

    total_memory_bytes = memory_for_n_bytes + memory_for_m_bytes + memory_for_output_bytes

    # The problem asks to output each number in the final equation.
    # The final answer format is x:o, so we print the values for x and o.
    # x = total_memory_bytes
    # o = max_covered_area
    print(f"{total_memory_bytes}:{max_covered_area}")

solve_tiling()