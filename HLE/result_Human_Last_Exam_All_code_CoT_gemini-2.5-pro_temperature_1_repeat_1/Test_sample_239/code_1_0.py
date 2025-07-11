def solve_theatre_square():
    """
    Calculates the number of flagstones needed for the given dimensions.
    """
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # To find the number of flagstones needed along one dimension, we need to
    # calculate the ceiling of the division. For positive integers, ceil(x/y)
    # can be calculated using integer division as (x + y - 1) // y.
    
    # Flagstones needed for the n-dimension
    flagstones_n = (n + a - 1) // a
    
    # Flagstones needed for the m-dimension
    flagstones_m = (m + a - 1) // a
    
    total_flagstones = flagstones_n * flagstones_m
    
    # The problem asks to output each number in the final equation.
    print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

solve_theatre_square()