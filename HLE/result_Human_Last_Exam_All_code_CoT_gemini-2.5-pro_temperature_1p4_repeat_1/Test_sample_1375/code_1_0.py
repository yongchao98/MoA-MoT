import math

def calculate_average_projection_area(k):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope,
    where n = 2k + 1.

    The formula used is: Area = (2 * sqrt(2k+1)) / (k!)^2
    """
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return

    # The dimension n is determined by k
    n = 2 * k + 1

    # Calculate the components of the formula
    try:
        k_factorial = math.factorial(k)
    except ValueError:
        print(f"Error: k={k} is too large to compute its factorial.")
        return
        
    k_factorial_sq = k_factorial ** 2
    sqrt_n = math.sqrt(n)
    
    # Assemble the final result
    numerator = 2 * sqrt_n
    denominator = k_factorial_sq
    result = numerator / denominator

    # Print the explanation and the calculation step-by-step
    print(f"For k = {k}:")
    print(f"The dimension of the cross-polytope is n = 2*k + 1 = {n}.")
    print("The formula for the average projection area is (2 * sqrt(n)) / (k!)^2.")
    print("\nCalculating the components:")
    print(f"sqrt(n) = sqrt({n}) \u2248 {sqrt_n:.4f}")
    print(f"k! = {k}! = {k_factorial}")
    print(f"(k!)^2 = ({k_factorial})^2 = {k_factorial_sq}")

    # Display the final equation with the computed numbers
    print("\nPlugging the values into the formula:")
    print(f"Average Area = (2 * {sqrt_n:.4f}) / {k_factorial_sq}")
    print(f"Final Result \u2248 {result}")

# --- Example Usage ---
# You can change the value of k here.
# For example, for a 3-dimensional cross-polytope (octahedron), n=3 so k=1.
# For a 5-dimensional cross-polytope, n=5 so k=2.
k_value = 3
calculate_average_projection_area(k_value)