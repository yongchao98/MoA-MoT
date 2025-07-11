import math

def solve_theatre_square():
    """
    Calculates the number of flagstones needed for the case where n = m = 4*10^9 and a = 1.
    This also demonstrates the logic for solving the problem without branching.
    """
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # The number of flagstones needed along one dimension is the ceiling of the length divided by the flagstone size.
    # For positive integers, ceil(x / y) can be calculated using integer division as (x + y - 1) // y.
    # This avoids floating-point inaccuracies and conditional branching.
    
    # Flagstones needed for the 'n' dimension
    flagstones_n = (n + a - 1) // a
    
    # Flagstones needed for the 'm' dimension
    flagstones_m = (m + a - 1) // a
    
    # Total flagstones is the product of the two dimensions.
    total_flagstones = flagstones_n * flagstones_m
    
    # Print each number in the final equation as requested.
    print(f"{int(flagstones_n)} * {int(flagstones_m)} = {total_flagstones}")

solve_theatre_square()
