def solve_tiling_problem():
    """
    This function calculates the maximal area of an NxM rectangle that can be
    covered with 2x1 bricks and the memory required by an equivalent C program.
    """
    
    # Given values for N and M
    n = 1000000001
    m = 1000000001

    # --- Calculate the maximal covered area (o) ---
    # The area of each brick is 2. The maximal covered area must be an even number.
    # We can find this by calculating the total area (n * m) and then finding the
    # largest even number less than or equal to it.
    # Using integer division, this is equivalent to (n * m) // 2 * 2.
    total_area = n * m
    max_covered_area = (total_area // 2) * 2

    # --- Calculate memory usage in a C program (x) ---
    # The values of n, m, and the output (~10^18) require a 64-bit integer type
    # in C, such as 'long long', which uses 8 bytes on most systems.
    # The problem asks for the memory used by variables n, m, and the output.
    memory_for_n = 8  # bytes for a 'long long'
    memory_for_m = 8  # bytes for a 'long long'
    memory_for_output = 8  # bytes for a 'long long'
    memory_used_bytes = memory_for_n + memory_for_m + memory_for_output

    # Print the result in the specified format 'x:o'.
    # The "final equation" is the output format itself, and this prints
    # both numbers required by it.
    print(f"{memory_used_bytes}:{max_covered_area}")

solve_tiling_problem()