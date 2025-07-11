import math

def solve():
    """
    This function calculates the least number of flagstones needed to pave a
    rectangular Theatre Square of size n*m with square flagstones of size a*a.
    """
    try:
        # Read the three integer inputs from one line, separated by spaces.
        n, m, a = map(int, input().split())

        # The number of flagstones needed along the n-dimension is ceil(n/a).
        # In integer arithmetic, this can be calculated as (n + a - 1) // a.
        # This formula avoids floating-point issues and branching (if statements).
        flagstones_n = (n + a - 1) // a
        
        # The number of flagstones needed along the m-dimension is ceil(m/a).
        flagstones_m = (m + a - 1) // a

        # The total number of flagstones is the product of the two dimensions.
        total_flagstones = flagstones_n * flagstones_m
        
        # Print the final equation as requested.
        print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

    except (ValueError, IndexError):
        print("Invalid input. Please provide three positive integers separated by spaces.")

# Execute the function to solve the problem.
# Example input: 6 7 4
# Expected output: 2 * 2 = 4
# Example input: 4000000000 4000000000 1
# Expected output: 4000000000 * 4000000000 = 16000000000000000000
solve()
