import math

def solve_alpha():
    """
    Calculates the value of alpha based on the problem derivation.
    """
    
    # Given constants
    T = math.log(10)
    B = 0.5 * (10**20) / (99**2)

    # Derived relationship between alpha and B
    # alpha = (sqrt(8*B) * (1 - e^(-3*T))) / 3

    # Step 1: Calculate sqrt(8*B)
    term1 = 8 * B
    sqrt_8B = math.sqrt(term1)
    
    # Step 2: Calculate 1 - e^(-3*T)
    # e^(-3*T) = e^(-3*ln(10)) = (e^ln(10))^(-3) = 10^(-3)
    exp_m3T = 10**(-3)
    term2 = 1 - exp_m3T

    # Step 3: Calculate the denominator
    denominator = 3.0

    # Final calculation for alpha
    alpha = (sqrt_8B * term2) / denominator
    
    # Print the equation with the calculated numbers
    print(f"The calculation for alpha is based on the formula: alpha = (sqrt(8*B) * (1 - e^(-3T))) / 3")
    print(f"Given T = ln(10), B = {B}")
    print(f"alpha = (sqrt(8 * {B}) * (1 - {exp_m3T})) / {denominator}")
    print(f"alpha = (sqrt({term1}) * ({term2})) / {denominator}")
    print(f"alpha = ({sqrt_8B} * {term2}) / {denominator}")
    print(f"alpha = {alpha}")

solve_alpha()