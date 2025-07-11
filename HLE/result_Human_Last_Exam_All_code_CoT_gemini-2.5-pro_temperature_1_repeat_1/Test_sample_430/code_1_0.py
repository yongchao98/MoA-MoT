import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site using a mean-field model
    based on a set of assumptions to resolve ambiguities in the problem statement.
    """
    # Parameters from the problem (in units of k_B*T)
    beta_eps1 = 0.1
    z_l = 4.0
    z_inter = 4.0

    # --- Assumptions to make the problem solvable ---
    # 1. Assume mu is negative to prevent infinite layer growth.
    beta_mu = -0.15
    # 2. Assume maximum number of layers K=2.
    K_max = 2
    # 3. Calculate lateral interaction energy eps_l based on K_max=2.
    beta_eps_l = (0.02)**K_max
    # 4. Assume eps_2 is related to eps_1 and z_inter.
    beta_eps2 = beta_eps1 / z_inter

    print("--- Model Parameters (in units of k_B*T) ---")
    print(f"epsilon_1 = {beta_eps1}")
    print(f"mu = {beta_mu} (Assumed negative sign)")
    print(f"epsilon_l = {beta_eps_l} (Assumed K_max={K_max})")
    print(f"epsilon_2 = {beta_eps2} (Assumed eps_2 = eps_1 / z_inter)")
    print(f"z_l = {z_l}")
    print(f"z_inter = {z_inter}")
    print("-------------------------------------------\n")

    # Iterative solver for self-consistent equations
    # Initial guess for layer coverages
    theta1, theta2 = 0.0, 0.0

    # Iterate to find the solution
    num_iterations = 10
    for i in range(num_iterations):
        # Calculate exponents for probability ratios
        arg1 = beta_eps1 + beta_mu + z_l * beta_eps_l * theta1
        arg2 = beta_eps2 + beta_mu + z_l * beta_eps_l * theta2
        
        # Calculate probability ratios x_i = p_i / p_{i-1}
        x1 = math.exp(arg1)
        x2 = math.exp(arg2)

        # Solve for probabilities p0, p1, p2 using the normalization condition
        # p0 * (1 + x1 + x1*x2) = 1
        p0 = 1.0 / (1.0 + x1 + x1 * x2)
        p1 = p0 * x1
        p2 = p0 * x1 * x2

        # Update layer coverages
        theta1 = p1 + p2
        theta2 = p2

    # The final result is the sum of the converged layer coverages
    avg_k = theta1 + theta2
    
    print("--- Converged Results ---")
    print(f"Layer 1 coverage (theta_1): {theta1:.4f}")
    print(f"Layer 2 coverage (theta_2): {theta2:.4f}\n")

    print("--- Final Calculation ---")
    print("The average number of adsorbed layers per site is <k> = theta_1 + theta_2")
    print(f"<k> = {theta1:.4f} + {theta2:.4f}")
    print(f"<k> = {avg_k:.4f}")

    return avg_k

# Run the calculation and store the final answer
final_answer = solve_adsorption()

# Final answer in the required format
# print(f"\n<<<{final_answer:.4f}>>>")