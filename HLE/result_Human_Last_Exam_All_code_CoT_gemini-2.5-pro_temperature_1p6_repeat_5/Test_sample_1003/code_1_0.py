import math

def solve_star_angle_ratio():
    """
    This function calculates the specified ratio of star angles based on the principles of special relativity.
    It follows the logical steps outlined above to find the ratio of Doppler factors,
    which is equivalent to the requested ratio of angle expressions.
    """

    # The problem reduces to finding the ratio D3/D1, which we found to be 1 / (1 + 1/sqrt(2)).
    
    # Define the numbers in the equation for the ratio R = n1 / (d1 + n2 / sqrt(d2))
    n1 = 1
    d1 = 1
    n2 = 1
    d2 = 2
    
    # Calculate the value of the expression
    # R = 1 / (1 + 1/sqrt(2))
    ratio_value = n1 / (d1 + n2 / math.sqrt(d2))
    
    # The simplified exact form is 2 - sqrt(2)
    simplified_value = 2 - math.sqrt(d2)

    # Print the explanation and the result.
    print("Based on the principles of special relativity, the expression (1 - cos(theta_14)) / (1 - cos(theta_34))")
    print("simplifies to the ratio of Doppler factors D_3 / D_1.")
    print("\nThis ratio is derived from the given angles in the second reference frame.")
    print(f"The equation for the ratio is: {n1} / ({d1} + {n2}/math.sqrt({d2}))")
    print(f"\nThis expression simplifies to the exact value of 2 - sqrt(2).")
    print(f"The numerical value is approximately: {simplified_value}")

solve_star_angle_ratio()
