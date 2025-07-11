def solve_tiling_problem():
    """
    Calculates the maximal covered area of an NxM rectangle with 2x1 bricks
    and the memory usage of an equivalent optimal C program.
    """
    
    # The given dimensions for the rectangle.
    n = 1000000001
    m = 1000000001

    # Step 1: Calculate the memory usage 'x' for an optimal C program.
    # n (long long) = 8 bytes
    # m (long long) = 8 bytes
    # output (__int128_t) = 16 bytes (since n*m overflows a 64-bit integer)
    # Total memory = 8 + 8 + 16 = 32 bytes.
    memory_usage_bytes = 32

    # Step 2: Calculate the maximal covered area 'o'.
    # Since both n and m are odd, the product n*m is odd.
    # The maximal covered area is the largest even number <= n*m, which is n*m - 1.
    # Python's integers handle arbitrary size, so overflow is not an issue here.
    max_area = n * m - 1
    
    # The instruction "Remember in the final code you still need to output each number 
    # in the final equation!" is interpreted as printing the final answer 'o'
    # which is the result of the calculation, along with 'x'.
    # The final output format is specified as x:o.
    print(f"{memory_usage_bytes}:{max_area}")

solve_tiling_problem()