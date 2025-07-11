import math

def calculate_average_projection_area(n):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    onto a fixed hyperplane, for odd n.

    Args:
        n: The dimension, an odd integer >= 3.
    """
    if not isinstance(n, int) or n % 2 == 0 or n < 3:
        print("Error: The dimension 'n' must be an odd integer greater than or equal to 3.")
        return

    # For n = 2k + 1, we find k
    k = (n - 1) // 2

    # The simplified formula for the average area is sqrt(n) / (k! * (k-1)!)
    # also written as sqrt(2k+1) / (k! * (k-1)!)

    # Calculate the components of the formula
    numerator = math.sqrt(n)
    k_factorial = math.factorial(k)
    k_minus_1_factorial = math.factorial(k - 1)
    denominator = k_factorial * k_minus_1_factorial

    # Calculate the final average area
    average_area = numerator / denominator

    # Output the steps of the calculation as requested
    print(f"For the {n}-dimensional cross-polytope:")
    print(f"We have n = 2k + 1, which gives k = (n-1)/2 = ({n}-1)/2 = {k}.")
    print("\nThe formula for the average projection area is: sqrt(n) / (k! * (k-1)!)")
    print("\nCalculation steps:")
    print(f"Average Area = sqrt({n}) / ({k}! * ({k}-1)!)")
    print(f"             = {numerator:.4f} / ({k_factorial} * {k_minus_1_factorial})")
    print(f"             = {numerator:.4f} / {denominator}")
    print(f"             = {average_area:.4f}")

# --- Example Usage ---
# You can change the value of n here to any odd integer >= 3.
# For example, n=3 is the octahedron. n=5, n=7 etc.
n_dimension = 5
calculate_average_projection_area(n_dimension)