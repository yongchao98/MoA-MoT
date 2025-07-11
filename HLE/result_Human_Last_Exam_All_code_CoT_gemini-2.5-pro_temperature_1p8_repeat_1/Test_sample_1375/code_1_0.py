import math

def average_projection_area(k):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope
    onto a fixed hyperplane, after a uniformly random rotation.
    The dimension n is given by n = 2k+1, where k is an integer >= 1.
    """
    if not isinstance(k, int) or k < 1:
        print("Error: k must be an integer greater than or equal to 1.")
        return

    n = 2 * k + 1
    
    # Calculate numerator
    numerator_val = 2 * math.sqrt(n)
    
    # Calculate denominator
    k_factorial = math.factorial(k)
    k_minus_1_factorial = math.factorial(k - 1)
    denominator_val = k_factorial * k_minus_1_factorial
    
    # Calculate the average area
    average_area = numerator_val / denominator_val
    
    # Output the result step-by-step
    print(f"For k = {k}, the dimension is n = 2*k + 1 = {n}.")
    print("The formula for the average projected area is: 2 * sqrt(2*k + 1) / (k! * (k-1)!)")
    print("\nCalculation:")
    print(f"Numerator = 2 * sqrt({n}) = {numerator_val:.4f}")
    print(f"Denominator = {k}! * {k-1}! = {k_factorial} * {k_minus_1_factorial} = {denominator_val}")
    print(f"\nAverage Area = {numerator_val:.4f} / {denominator_val} = {average_area:.4f}")

# Example usage for k=3 (which corresponds to n=7)
average_projection_area(3)
