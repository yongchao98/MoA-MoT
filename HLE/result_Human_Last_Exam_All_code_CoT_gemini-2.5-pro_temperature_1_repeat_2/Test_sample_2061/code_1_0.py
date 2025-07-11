import math

def solve_for_alpha():
    """
    This function calculates the value of alpha based on the derived formula.
    """
    # Given values from the problem statement
    # T = ln(10)
    # B = (1/2) * 10^20 / 99^2
    T_val = math.log(10)
    B_val = (1/2) * (10**20) / (99**2)

    # Derived formula for alpha:
    # alpha = (sqrt(8 * B) * (1 - exp(-3 * T))) / 3
    
    # Constants in the formula
    c1 = 8
    c2 = 1
    c3 = 3
    c4 = 3

    # Step-by-step calculation
    sqrt_term = math.sqrt(c1 * B_val)
    exp_term = math.exp(-c3 * T_val)
    numerator = sqrt_term * (c2 - exp_term)
    alpha = numerator / c4
    
    # Output the details of the calculation
    print("Based on the problem analysis, the formula for alpha is:")
    print(f"alpha = (sqrt(c1 * B) * (c2 - exp(-c3 * T))) / c4")
    print("\nSubstituting the given values:")
    print(f"c1 = {c1}")
    print(f"B = {B_val}")
    print(f"c2 = {c2}")
    print(f"c3 = {c3}")
    print(f"T = {T_val}")
    print(f"c4 = {c4}")
    
    print("\nResult of the calculation:")
    print(f"alpha = {alpha}")

solve_for_alpha()