def solve_tiling_problem():
    """
    Calculates the maximal covered area of a rectangle with 2x1 bricks and
    the memory usage of an equivalent optimal C program.
    """
    
    # Input dimensions as per the problem statement.
    # We use Python's arbitrary-precision integers, which handle these large numbers automatically.
    N = 1000000001
    M = 1000000001
    
    # The maximal area is the total area of the rectangle, rounded down to the nearest even number.
    # This can be calculated using integer division.
    # The number of bricks is (N * M) // 2
    # The covered area is the number of bricks * 2
    total_area = N * M
    max_covered_area = (total_area // 2) * 2

    # In C, to store N, M, and the result (which is ~10^18), we would need 64-bit integers (`long long`).
    # On most systems, a `long long` is 8 bytes.
    # The problem asks for memory used by variables (n, m, the output).
    # Memory for n (long long): 8 bytes
    # Memory for m (long long): 8 bytes
    # Memory for the output (long long): 8 bytes
    # Total memory = 8 + 8 + 8 = 24 bytes.
    memory_in_bytes = 24
    
    # The final output format is x:o
    # x is the memory used (in bytes)
    # o is the maximal area
    print(f"{memory_in_bytes}:{max_covered_area}")

solve_tiling_problem()