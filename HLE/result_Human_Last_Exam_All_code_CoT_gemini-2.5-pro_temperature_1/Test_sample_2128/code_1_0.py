import math

def solve():
    """
    Calculates the value of 1/p_1000 based on the derived formula.
    """
    n = 1000
    
    # The formula for 1/p_n is 4 * cos^2(pi / (n + 2))
    # We need to calculate this for n = 1000.
    
    # The final equation is: 1/p_1000 = 4 * cos^2(pi / 1002)
    # The numbers in this equation are:
    constant_factor = 4
    numerator_of_angle = math.pi
    denominator_of_angle = n + 2

    print(f"To find 1/p_n for n = {n}, we use the formula: C * cos^2(Num / Den)")
    print("The numbers in the final equation are:")
    print(f"C (Constant Factor) = {constant_factor}")
    print(f"Num (Numerator of angle) = pi = {numerator_of_angle}")
    print(f"Den (Denominator of angle) = n + 2 = {denominator_of_angle}")

    # Perform the calculation
    angle = numerator_of_angle / denominator_of_angle
    cos_squared = math.cos(angle)**2
    result = constant_factor * cos_squared
    
    print(f"\nThe value of 1/p_{n} is: {result}")

solve()