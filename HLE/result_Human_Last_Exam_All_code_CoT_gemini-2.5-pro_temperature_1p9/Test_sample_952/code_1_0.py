import math

def solve_diamond_puzzle():
    """
    This function calculates the largest value k based on the grid size N.
    k represents the largest number of diamonds such that any arrangement of k diamonds
    is guaranteed to have at least one diamond that can be moved to an adjacent cell.
    
    The logic is k = M_min - 1, where M_min is the minimum number of diamonds
    in a "frozen" configuration. For an N x N grid, M_min = ceil(N^2 / 3).
    """
    N = 2024
    
    # Calculate N squared
    N_squared = N * N
    
    # M_min is the ceiling of N_squared / 3.
    # In integer arithmetic, ceil(a/b) can be calculated as (a + b - 1) // b.
    M_min = (N_squared + 3 - 1) // 3
    
    # The answer k is M_min - 1
    k = M_min - 1
    
    # The final equation is k = ceil(N * N / 3) - 1
    # We will print the numbers that form this equation.
    print(f"For a grid of size {N} x {N}:")
    print(f"The minimum number of diamonds in a frozen arrangement (M_min) is ceil({N} * {N} / 3) = {M_min}")
    print(f"The largest value k is M_min - 1.")
    print(f"Final calculation: {k} = {M_min} - 1")
    print(f"The value is {k}")

solve_diamond_puzzle()