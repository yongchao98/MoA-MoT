import math

def solve_alpha():
    """
    This function solves for the value of alpha based on the problem description.
    """
    # Given parameters
    T = math.log(10)
    B = (1/2) * (10**20) / (99**2)
    
    # Step 1: From the integral evaluation, we derived K^2 / 8 = B.
    # Solve for K. Since alpha > 0, K must be positive.
    K_squared = 8 * B
    K = math.sqrt(K_squared)
    
    # Step 2: From the solvability condition, K = 3*alpha / (1 - e^(-3T)).
    # Solve for alpha.
    # First, calculate the denominator 1 - e^(-3T).
    # e^(-3T) = (e^T)^(-3) = 10^(-3) = 1/1000.
    denominator = 1 - math.exp(-3 * T)
    
    # alpha = K * (1 - e^(-3T)) / 3
    alpha = K * denominator / 3
    
    # Print the values in the final equation K = 3*alpha / (1 - e^(-3T))
    print(f"The equation relating K and alpha is K = 3*alpha / (1 - e^(-3T))")
    print(f"From B, we found K = {K}")
    print(f"The denominator (1 - e^(-3T)) for T=ln(10) is {denominator}")
    print(f"Substituting K and the denominator, we get:")
    print(f"{K} = 3 * alpha / {denominator}")
    print(f"Solving for alpha gives: alpha = {alpha}")

solve_alpha()