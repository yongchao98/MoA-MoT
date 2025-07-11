import math

def solve_theatre_square():
    """
    Calculates the least number of flagstones to pave a rectangular square.
    """
    try:
        # Prompt the user for input
        print("Enter n, m, and a, separated by spaces:")
        n, m, a = map(int, input().split())

        # Ensure inputs are positive as per the problem description
        if not (n > 0 and m > 0 and a > 0):
            print("Error: Inputs n, m, and a must be positive integers.")
            return

        # Calculate the number of flagstones needed for the n-side.
        # We use integer division to find the ceiling of a division:
        # ceil(x / y) is equivalent to (x + y - 1) // y in integer arithmetic.
        # This avoids floating-point inaccuracies with large numbers.
        stones_n = (n + a - 1) // a

        # Calculate the number of flagstones for the m-side.
        stones_m = (m + a - 1) // a

        # The total number of flagstones is the product of the two sides.
        total_stones = stones_n * stones_m

        # Print the final equation with all numbers, as requested.
        print(f"\nTo cover a {n}x{m} square with {a}x{a} flagstones:")
        print(f"Number of stones along the n-side: {stones_n}")
        print(f"Number of stones along the m-side: {stones_m}")
        print(f"Final Equation: {stones_n} * {stones_m} = {total_stones}")

    except ValueError:
        print("Invalid input. Please enter three integer numbers separated by spaces.")
    except Exception as e:
        print(f"An error occurred: {e}")

# Execute the function
solve_theatre_square()