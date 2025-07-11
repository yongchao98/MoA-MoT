import math

def solve_theatre_square():
    """
    This function solves the Theatre Square problem for given inputs and
    prints the final equation.
    """
    # The input values for n, m, and a.
    # Using the large values from Question 4 as an example.
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # Calculate the number of flagstones needed along the width (n).
    # This is equivalent to ceil(n/a). We can use integer arithmetic
    # to avoid floating-point numbers: (n + a - 1) // a
    flagstones_n = (n + a - 1) // a
    
    # Calculate the number of flagstones needed along the length (m).
    # This is equivalent to ceil(m/a).
    flagstones_m = (m + a - 1) // a
    
    # The total number of flagstones is the product of the two dimensions.
    total_flagstones = flagstones_n * flagstones_m

    # Per the instructions, output each number in the final equation.
    # The equation represents (flagstones for n-side) * (flagstones for m-side) = total flagstones.
    print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

solve_theatre_square()