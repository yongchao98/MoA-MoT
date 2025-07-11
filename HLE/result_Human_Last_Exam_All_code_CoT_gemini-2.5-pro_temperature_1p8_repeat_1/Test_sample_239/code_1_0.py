import math

def solve_theatre_square():
    """
    Calculates the minimum number of flagstones to pave the Theatre Square
    for the specific inputs given in the problem.
    """
    # Input values from question 4
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # Using math.ceil to find the number of flagstones needed for each dimension.
    # This is equivalent to the (x + y - 1) // y integer math trick and avoids 'if' statements.
    # We cast to int because the inputs are whole numbers and the result of ceil is a float.
    na = int(math.ceil(n / a))
    ma = int(math.ceil(m / a))

    # Calculate the total number of flagstones
    total_flagstones = na * ma

    # Print the final equation as requested
    print(f"{na} * {ma} = {total_flagstones}")

solve_theatre_square()