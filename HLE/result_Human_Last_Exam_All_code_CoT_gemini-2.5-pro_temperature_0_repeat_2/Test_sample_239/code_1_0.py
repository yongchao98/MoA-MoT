import math

def solve_theatre_square():
    """
    Calculates the minimum number of flagstones to pave a rectangular Theatre Square.
    """
    try:
        # 1. Read the three positive integers n, m, and a from a single line of input.
        n, m, a = map(int, input("Enter n, m, and a (space-separated): ").split())

        # 2. Calculate the number of flagstones needed for the n-dimension.
        # This is the ceiling of n/a. Using integer arithmetic (n + a - 1) // a
        # is robust and avoids floating-point inaccuracies.
        flagstones_n = (n + a - 1) // a

        # 3. Calculate the number of flagstones needed for the m-dimension.
        flagstones_m = (m + a - 1) // a

        # 4. The total number of flagstones is the product of the two dimensions.
        total_flagstones = flagstones_n * flagstones_m

        # 5. Print the final equation showing each number, as requested.
        print(f"\nTo cover a {n}x{m} square with {a}x{a} flagstones, the calculation is:")
        print(f"({n} + {a} - 1) // {a}  *  ({m} + {a} - 1) // {a}")
        print(f"= {flagstones_n} * {flagstones_m}")
        print(f"= {total_flagstones}")

    except (ValueError, IndexError):
        print("Invalid input. Please enter three space-separated integer numbers.")
    except ZeroDivisionError:
        print("Error: The flagstone size 'a' cannot be zero.")

# Execute the function
# You can test this with the values from the problem, e.g., "4000000000 4000000000 1"
solve_theatre_square()