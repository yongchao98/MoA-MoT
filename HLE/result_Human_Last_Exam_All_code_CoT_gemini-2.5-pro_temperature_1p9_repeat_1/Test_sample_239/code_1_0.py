import math

def solve_flagstones():
    """
    Calculates the minimum number of flagstones to cover a rectangular area.
    This function implements the branchless solution discussed in the analysis.
    """
    try:
        # Read the three integer numbers from a single line of input
        n, m, a = map(int, input().split())

        # Validate inputs based on the problem description
        if not (0 < n <= 4 * 10**9 and 0 < m <= 4 * 10**9 and 0 < a < 10**6):
            print("Input values are out of the specified range.")
            return

        # Calculate the number of flagstones needed along the n and m dimensions
        # using ceiling division. In integer arithmetic, ceil(x/y) is (x + y - 1) // y.
        # Python's math.ceil returns a float, so we use integer math for precision.
        # Python integers support arbitrary size, so overflow is not an issue.
        
        num_a_for_n = (n + a - 1) // a
        num_a_for_m = (m + a - 1) // a

        # Total number of flagstones is the product of the two dimensions
        total_flagstones = num_a_for_n * num_a_for_m

        # Output the numbers in the final equation as requested
        print(f"{num_a_for_n} * {num_a_for_m} = {total_flagstones}")

    except (ValueError, IndexError):
        print("Invalid input. Please provide three space-separated integers (n m a).")

# To run this code, you would execute the script and then type the inputs,
# for example: 4000000000 4000000000 1
if __name__ == '__main__':
    solve_flagstones()
