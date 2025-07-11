import math

def calculate_average_projection_area(n):
    """
    Calculates the average projection area of an n-dimensional cross-polytope,
    where n is an odd integer.
    """
    if not isinstance(n, int) or n < 1 or n % 2 == 0:
        print("Error: The dimension n must be a positive odd integer.")
        return

    if n == 1:
        print("For n = 1, the cross-polytope is the interval [-1, 1].")
        print("Its projection onto a 0-dimensional hyperplane (a point) is a point.")
        print("The 0-dimensional area is 0.")
        return

    k = (n - 1) // 2
    
    # Calculate components of the formula: E_n = (2 * sqrt(n)) / (n * k! * (k-1)!)
    sqrt_n = math.sqrt(n)
    k_fact = math.factorial(k)
    k_minus_1_fact = math.factorial(k - 1)
    
    numerator = 2 * sqrt_n
    denominator = n * k_fact * k_minus_1_fact
    
    result = numerator / denominator

    print(f"The average area of a projection of the {n}-dimensional cross-polytope is calculated.")
    print(f"Given dimension n = {n}.")
    print(f"We have k = (n-1)/2 = {k}.")
    print("\nThe formula for the average area is:")
    print("  E_n = (2 * sqrt(n)) / (n * k! * (k-1)!)")
    
    print("\nSubstituting the values:")
    print(f"  E_n = (2 * sqrt({n})) / ({n} * {k}! * {k-1}!)")
    print(f"  E_n = (2 * {sqrt_n:.5f}) / ({n} * {k_fact} * {k_minus_1_fact})")
    print(f"  E_n = {numerator:.5f} / {denominator:.5f}")
    print(f"  E_n = {result:.5f}")

# Example calculation for n=5
n_dimension = 5
calculate_average_projection_area(n_dimension)
