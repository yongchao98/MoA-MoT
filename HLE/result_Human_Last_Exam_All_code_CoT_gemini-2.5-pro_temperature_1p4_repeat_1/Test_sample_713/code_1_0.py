import math

def calculate_polytope_volume_factor(d):
    """
    Calculates the factor relating the volume of the polytope P to the simplex T.

    Args:
        d (int): The dimension of the space. Must be a non-negative integer.

    Returns:
        None. Prints the relationship.
    """
    if not isinstance(d, int) or d < 1:
        print("Error: The dimension d must be an integer greater than or equal to 1.")
        return

    # For a d-dimensional simplex T with volume V, the volume of the
    # associated polytope P is given by the formula: Vol(P) = d! * V.
    # We calculate d! (d factorial).
    try:
        factor = math.factorial(d)
    except ValueError:
        print(f"Error: Factorial is not defined for negative dimension {d}.")
        return


    print(f"For a non-degenerate simplex T in {d} dimensions with volume V,")
    print(f"the volume of the corresponding polytope P is always {factor} times the volume of T.")
    print("\nThis gives the equation:")
    # We output each number in the final equation as requested.
    # Let's represent the variables symbolically in the printout.
    print(f"Vol(P) = {d}! * Vol(T) = {factor} * Vol(T)")

if __name__ == '__main__':
    try:
        # You can change the dimension d here.
        d_dimension = 3
        # For a more interactive experience, you could use:
        # d_dimension = int(input("Enter the dimension d: "))
        calculate_polytope_volume_factor(d_dimension)
    except (ValueError, TypeError):
        print("Invalid input. Please enter an integer for the dimension.")
