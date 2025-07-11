import math

def solve_theatre_square():
    """
    Calculates the minimum number of flagstones to pave a rectangular square.
    """
    try:
        # Read n, m, and a from a single line of input
        n, m, a = map(int, input().split())

        # Check for valid positive integer inputs as per the problem description
        if not (n > 0 and m > 0 and a > 0):
            print("Inputs n, m, and a must be positive integers.")
            return

        # Calculate the number of flagstones needed for the length (n)
        # This is equivalent to ceil(n/a)
        flagstones_n = (n + a - 1) // a

        # Calculate the number of flagstones needed for the width (m)
        # This is equivalent to ceil(m/a)
        flagstones_m = (m + a - 1) // a

        # Calculate the total number of flagstones
        total_flagstones = flagstones_n * flagstones_m

        # Print the final equation with the computed values
        # This fulfills the requirement: "output each number in the final equation"
        print(f"{flagstones_n} * {flagstones_m} = {total_flagstones}")

    except (ValueError, IndexError):
        print("Invalid input. Please provide three space-separated integers (n m a).")

# Execute the function
solve_theatre_square()