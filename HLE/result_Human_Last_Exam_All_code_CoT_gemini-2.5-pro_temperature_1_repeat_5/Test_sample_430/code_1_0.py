import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers based on the
    interpreted model and self-consistency argument.
    """
    # Parameters from the problem statement, normalized by k_B*T
    mu_norm = 0.15
    eps1_norm = 0.1
    z_l = 4
    
    # Based on the self-consistency argument, the problem reduces to k=1
    k = 1
    
    # Calculate epsilon_l for k=1
    epsl_norm = (0.02)**k
    
    # The problem reduces to solving the Fowler-Guggenheim equation for theta:
    # ln(theta / (1 - theta)) = mu_norm + eps1_norm + z_l * epsl_norm * theta
    # Let C1 = mu_norm + eps1_norm
    # Let C2 = z_l * epsl_norm
    C1 = mu_norm + eps1_norm
    C2 = z_l * epsl_norm

    print("--- Model Parameters and Equation ---")
    print(f"Based on a self-consistency argument, we solve for the k=1 case.")
    print(f"Maximum layers (k) = {k}")
    print(f"Chemical Potential (mu / k_B*T) = {mu_norm}")
    print(f"First Layer Energy (epsilon_1 / k_B*T) = {eps1_norm}")
    print(f"Lateral Coordination Number (z_l) = {z_l}")
    print(f"Lateral Interaction Energy (epsilon_l / k_B*T) = (0.02)^k = {epsl_norm:.2f}")
    
    print("\n--- Equation to Solve ---")
    print("ln(theta / (1 - theta)) = (mu + epsilon_1)/(k_B*T) + (z_l * epsilon_l / (k_B*T)) * theta")
    print("Which simplifies to:")
    print(f"ln(theta / (1 - theta)) = {C1:.2f} + {C2:.2f} * theta")

    # Solve the transcendental equation iteratively for theta
    # theta = 1 / (1 + exp(-(C1 + C2 * theta)))
    theta = 0.5  # Initial guess
    for _ in range(20): # 20 iterations is more than enough for convergence
        theta = 1.0 / (1.0 + math.exp(-(C1 + C2 * theta)))
        
    print("\n--- Result ---")
    print(f"The average number of adsorbed layers per site, <k>, is {theta:.4f}")
    
    # Final answer in the required format
    print(f"\n<<<{theta:.4f}>>>")

solve_adsorption()