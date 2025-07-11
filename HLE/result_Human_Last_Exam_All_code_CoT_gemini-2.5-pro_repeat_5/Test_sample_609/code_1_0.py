import math

def calculate_area_ratio(n):
    """
    Calculates the area ratio of an n-sided polygon constructed from a 2n-sided polygon.

    Args:
        n (int): The number of sides of the outer polygon (must be >= 3).
    """
    if n < 3:
        print("The number of sides 'n' must be 3 or greater.")
        return

    print(f"Calculating for a {n}-sided polygon constructed from a {2*n}-sided polygon.")
    print("-" * 50)
    
    # The general formula for the area ratio
    print("The general formula for the ratio is: tan(π/n) / (2 * tan(π/(2n)))")
    
    # Substitute the value of n into the formula description
    print(f"For n = {n}, the equation is: tan(π/{n}) / (2 * tan(π/({2*n})))")
    
    # Calculate the values
    pi_val = math.pi
    numerator_angle = pi_val / n
    denominator_angle = pi_val / (2 * n)
    
    numerator = math.tan(numerator_angle)
    denominator = 2 * math.tan(denominator_angle)
    
    ratio = numerator / denominator
    
    # Print the components of the final equation
    print("\nBreaking down the calculation:")
    print(f"  Numerator value: tan(π/{n}) = tan({round(numerator_angle, 4)}) = {numerator}")
    print(f"  Denominator value: 2 * tan(π/{2*n}) = 2 * tan({round(denominator_angle, 4)}) = {denominator}")
    
    # Print the final result
    print("\nFinal Result:")
    print(f"The area of the {n}-sided polygon is {ratio} times larger than the {2*n}-sided polygon.")
    
# Example from the problem description: n=3 (equilateral triangle from a hexagon)
calculate_area_ratio(3)

<<<1 / (2 * (1/cos(pi/n) - 1))>>>