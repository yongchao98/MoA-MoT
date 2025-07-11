import math

# Computer:
def solve_theatre_square():
    """
    Calculates the least number of flagstones to pave a rectangular square.
    This script demonstrates the correct logic using Python's arbitrary-precision integers,
    which do not suffer from the overflow issues of the XVM C `long` type.
    """
    # Inputs from the problem's fourth question
    n = 4 * 10**9
    m = 4 * 10**9
    a = 1

    # To find the number of flagstones for a dimension, we divide its length by the
    # flagstone side and round up. This is the ceiling function.
    # For positive integers x and y, ceil(x / y) can be calculated using integer
    # division as (x + a - 1) // a.

    # Calculate the number of flagstones needed for each dimension
    stones_n = (n + a - 1) // a
    stones_m = (m + a - 1) // a

    # The total number is the product of the two dimensions
    total_stones = stones_n * stones_m

    # Print the full equation and the final result as requested
    print(f"To cover a {n}x{m} area with {a}x{a} flagstones, the calculation is:")
    print(f"(({n} + {a} - 1) // {a}) * (({m} + {a} - 1) // {a}) = {stones_n} * {stones_m} = {total_stones}")

solve_theatre_square()