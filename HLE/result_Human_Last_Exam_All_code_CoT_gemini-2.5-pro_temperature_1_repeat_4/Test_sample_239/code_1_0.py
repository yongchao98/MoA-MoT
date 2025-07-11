import sys
import math

def solve_theatre_square():
    """
    This script solves the Theatre Square problem.

    It reads three space-separated integers (n, m, a) from a single line
    of standard input and calculates the minimum number of a*a flagstones
    to pave an n*m square.
    """
    try:
        # Read n, m, and a from a single line of input
        n, m, a = map(int, sys.stdin.readline().split())

        # The number of flagstones needed along one dimension is the ceiling
        # of the dimension's length divided by the flagstone's side length.
        # In integer arithmetic, ceil(x / y) can be calculated as (x + y - 1) // y.
        # This avoids floating point math and 'if' statements.
        
        # Calculate flagstones needed for the 'n' side
        tiles_n = (n + a - 1) // a
        
        # Calculate flagstones needed for the 'm' side
        tiles_m = (m + a - 1) // a
        
        # Total flagstones is the product of the two sides
        total_tiles = tiles_n * tiles_m

        # As requested, print the equation with each number included
        print(f"To cover a {n}x{m} square with {a}x{a} flagstones:")
        print(f"Tiles needed for side n: ({n} + {a} - 1) // {a} = {tiles_n}")
        print(f"Tiles needed for side m: ({m} + {a} - 1) // {a} = {tiles_m}")
        print(f"Total tiles required: {tiles_n} * {tiles_m} = {total_tiles}")

    except (ValueError, IndexError):
        print("Please provide three positive integers (n, m, a) separated by spaces on a single line.")
        print("Example: 6 7 4")

if __name__ == "__main__":
    solve_theatre_square()
