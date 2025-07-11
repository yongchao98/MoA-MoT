import math

def solve_theatre_square():
    """
    Calculates the number of flagstones for the case where n=m=4*10^9 and a=1.
    """
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # To find the number of flagstones needed along one dimension, we calculate
    # the ceiling of the division of the dimension length by the flagstone side length.
    # In integer arithmetic, ceil(x / y) is equivalent to (x + y - 1) // y.
    flagstones_n = (n + a - 1) // a
    flagstones_m = (m + a - 1) // a

    # The total number of flagstones is the product of the number needed for each dimension.
    total_flagstones = flagstones_n * flagstones_m

    # Per the instruction to "output each number in the final equation", we format the print statement this way.
    print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

solve_theatre_square()