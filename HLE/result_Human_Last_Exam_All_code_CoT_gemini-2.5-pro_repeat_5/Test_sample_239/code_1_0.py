import sys

# This script solves the "Theatre Square" problem by calculating the least number of flagstones.
# It reads three positive integer numbers n, m, and a from a single line of standard input.

def solve():
    """
    Reads input, calculates the number of flagstones, and prints the result.
    """
    try:
        # Read n, m, and a from a single line of input.
        n, m, a = map(int, sys.stdin.readline().split())
    except (ValueError, IndexError):
        # Handle cases with incorrect input format.
        print("Invalid input. Please provide three positive integers n, m, a separated by spaces on a single line.", file=sys.stderr)
        return

    # To pave an area of size 'n' with stones of size 'a', we need ceil(n/a) stones.
    # In integer arithmetic, the ceiling of a division (ceil(x/y)) can be calculated
    # efficiently without branching or floating-point math using the formula: (x + y - 1) // y.

    # Calculate the number of flagstones needed for the n-dimension (width).
    flagstones_n = (n + a - 1) // a

    # Calculate the number of flagstones needed for the m-dimension (length).
    flagstones_m = (m + a - 1) // a

    # The total number of flagstones is the product of the quantities for each dimension.
    total_flagstones = flagstones_n * flagstones_m

    # Print the final equation as requested, showing each number involved.
    print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

if __name__ == "__main__":
    solve()