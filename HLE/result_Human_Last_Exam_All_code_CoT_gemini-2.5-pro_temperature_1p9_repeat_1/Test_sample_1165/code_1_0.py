import math

def solve_ode_fluctuations():
    """
    This function calculates and prints the scaling relation for R(epsilon)
    for both the uniform and normal distribution cases based on the analytical
    derivation outlined in the plan.
    """

    print("--- Analysis of fluctuations for the given ODE ---")

    # Part 1: z_i are ordered uniform random variables
    
    # According to the derivation, for the uniform case, the maximum variance is
    # max(Var[y-y0]) = epsilon^4 * max(Var[sum(G)])
    # where max(Var[sum(G)]) is approximately 1 / (48 * epsilon^3).
    # So, R^2 = epsilon^4 / (48 * epsilon^3) = epsilon / 48.
    
    # R = sqrt(epsilon / 48) = (1 / sqrt(48)) * epsilon^0.5
    # 1/sqrt(48) = 1 / (4 * sqrt(3))
    
    C1 = 1 / (4 * math.sqrt(3))
    power1 = 0.5
    
    print("\n[Case 1: z_i ~ Uniform([0, epsilon^-1])]")
    print("The estimated maximum magnitude of fluctuations is R = (max|Var[y(x) - y0(x)]|)^1/2.")
    print("The scaling relation is derived as:")
    print(f"R(epsilon) = (1 / (4 * sqrt(3))) * epsilon^{power1}")
    print("Numerically, this is:")
    print(f"R(epsilon) = {C1:.4f} * epsilon^{power1}")


    # Part 2: z_i are independent normal random variables
    
    # For the normal case, z_i ~ Normal(i, 0.5), so sigma = 0.5.
    # The derivation shows that R^2 is approximately
    # epsilon^4 * (sigma^2 / (4*epsilon)).
    # R^2 = (sigma^2 / 4) * epsilon^3
    # R = (sigma / 2) * epsilon^1.5

    sigma = 0.5
    C2 = sigma / 2
    power2 = 1.5
    
    print("\n[Case 2: z_i ~ Normal(i, 0.5)]")
    print("If the source locations are changed, the scaling relation also changes.")
    print("The new scaling relation is derived as:")
    print(f"R(epsilon) = ({sigma} / 2) * epsilon^{power2}")
    print("Numerically, this is:")
    print(f"R(epsilon) = {C2:.4f} * epsilon^{power2}")
    
    print("\nComparing the two results, the scaling for R(epsilon) is not expected to remain the same.")

solve_ode_fluctuations()