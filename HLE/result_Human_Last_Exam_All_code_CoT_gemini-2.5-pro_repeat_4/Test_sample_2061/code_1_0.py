import math

def solve_alpha():
    """
    This function solves for the value of alpha based on the problem description.
    It follows the plan outlined above, assuming a typo in the definition of A (A^4 instead of A).
    """
    # Step 1: Define the given constants
    T = math.log(10)
    B = 0.5 * (10**20) / (99**2)

    # Step 2: Calculate intermediate exponential terms
    e_m2T = math.exp(-2 * T)
    e_m3T = math.exp(-3 * T)

    # Step 3: Use the derived formula to calculate alpha^2
    # The formula is alpha^2 = B * 2 * (1 - e^{-3T}) * (1 - e^{-2T})^2 / (2 - 3e^{-2T} + e^{-3T})
    
    val_1_minus_e_m2T = 1 - e_m2T
    val_1_minus_e_m3T = 1 - e_m3T
    
    numerator = 2 * B * val_1_minus_e_m3T * (val_1_minus_e_m2T**2)
    denominator = 2 - 3 * e_m2T + e_m3T
    
    alpha_squared = numerator / denominator
    
    # Since alpha > 0
    alpha = math.sqrt(alpha_squared)
    
    print(f"The calculated value of alpha is: {alpha}\n")

    # Step 4: Verify the result by plugging alpha back into the integral equation
    # The equation is K * A^4 / 4 - A^8 / 8 = B
    
    # Calculate K and A^4 using our found alpha
    K = (3 * alpha) / val_1_minus_e_m3T
    A_pow_4 = (2 * alpha) / val_1_minus_e_m2T
    A_pow_8 = A_pow_4**2
    
    # Calculate the left-hand side of the equation
    LHS = K * A_pow_4 / 4 - A_pow_8 / 8

    # Output the numbers in the final equation
    print("Verification of the solution:")
    print("The final equation is: K * A^4 / 4 - A^8 / 8 = B")
    print(f"Using the calculated alpha = {alpha}:")
    print(f"K = {K}")
    print(f"A^4 = {A_pow_4}")
    print(f"A^8 = {A_pow_8}")
    print(f"LHS = ({K}) * ({A_pow_4}) / 4 - ({A_pow_8}) / 8 = {LHS}")
    print(f"B = {B}")
    print(f"Does LHS equal B? {math.isclose(LHS, B)}")

solve_alpha()