def solve_tiling_problem():
    """
    Calculates the maximal area of an NxM rectangle that can be covered by 2x1 bricks
    and analyzes the memory usage of an equivalent C program.
    """
    # Input values for the rectangle dimensions
    N = 1000000001
    M = 1000000001

    # The maximum number of 2x1 bricks is the total area integer-divided by 2.
    # The maximum covered area is the number of bricks multiplied by 2.
    # Python's integers handle arbitrary size, so there is no overflow risk here.
    maximal_area = (N * M) // 2 * 2

    # In a C program, N and M could fit in a 4-byte 'int'.
    # However, their product (N * M) is ~10^18, which requires a 64-bit integer type
    # ('long long', 8 bytes) to avoid overflow.
    # To perform the calculation safely, n, m, and the result variable would all
    # be declared as 'long long'.
    # Memory for n (long long): 8 bytes
    # Memory for m (long long): 8 bytes
    # Memory for output (long long): 8 bytes
    memory_in_bytes = 8 + 8 + 8

    # The final output should be in the format x:o
    # x = memory used in bytes
    # o = the maximal area
    print(f"{memory_in_bytes}:{maximal_area}")

solve_tiling_problem()