import sys

def solve_theatre_square():
    """
    Reads three integers n, m, and a from a single line of input,
    and calculates the minimum number of a*a flagstones to pave
    an n*m square.
    """
    try:
        # Read n, m, and a from a single line of space-separated input.
        n, m, a = map(int, sys.stdin.readline().split())

        # Calculate the number of flagstones needed along the length 'n'.
        # This is equivalent to the ceiling of n/a, which can be calculated
        # using integer division as (n + a - 1) // a.
        flagstones_n = (n + a - 1) // a

        # Calculate the number of flagstones needed along the width 'm'.
        flagstones_m = (m + a - 1) // a

        # The total number of flagstones is the product of the two dimensions.
        total_flagstones = flagstones_n * flagstones_m

        # Print the numbers used in the final equation as requested.
        print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

    except (ValueError, IndexError):
        print("Invalid input. Please provide three positive integers separated by spaces on a single line.")

# Execute the main function to solve the problem
solve_theatre_square()