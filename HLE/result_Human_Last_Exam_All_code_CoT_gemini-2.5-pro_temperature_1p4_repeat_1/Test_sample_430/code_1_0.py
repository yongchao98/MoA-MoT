import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site using a
    mean-field approximation for a two-layer system.
    """
    # 1. Define parameters from the problem statement
    # Note: All energies are given in units of k_B*T, so we can work with their
    # values multiplied by beta = 1/(k_B*T). The temperature T=318K is not explicitly needed.
    
    beta_mu = 0.15  # mu / (k_B*T)
    beta_eps1 = 0.1 # epsilon_1 / (k_B*T)
    zl = 4          # Lateral coordination number
    z_inter = 4     # Vertical coordination number
    
    # 2. Make assumptions based on the problem statement
    # Assume k_max = 2 (simplest non-trivial case for multiple layers)
    k_max = 2
    
    # Calculate beta*epsilon_l based on k_max=2
    beta_eps_l = 0.02**k_max
    
    # Calculate beta*epsilon_2 based on the assumption epsilon_2 = z_inter * epsilon_l
    beta_eps2 = z_inter * beta_eps_l
    
    # 3. Solve the self-consistent equations by iteration
    
    # Initial guess for layer coverages
    theta_1 = 0.0
    theta_2 = 0.0
    
    # Iteration loop for convergence
    for _ in range(20): # 20 iterations are more than enough for convergence
        # Energy terms (depend on theta_1 and theta_2)
        # beta*E_m = -beta*sum(eps_j) - beta*zl*eps_l*sum(theta_j)
        beta_E1 = -beta_eps1 - zl * beta_eps_l * theta_1
        beta_E2 = -(beta_eps1 + beta_eps2) - zl * beta_eps_l * (theta_1 + theta_2)

        # Statistical weights for states m=0, 1, 2
        W0 = 1.0
        W1 = math.exp(1 * beta_mu - beta_E1)
        W2 = math.exp(2 * beta_mu - beta_E2)

        # Partition function for a single site
        Z_site = W0 + W1 + W2

        # Update theta values based on the new probabilities
        # theta_1 = p1 + p2
        # theta_2 = p2
        theta_1_new = (W1 + W2) / Z_site
        theta_2_new = W2 / Z_site
        
        # Check for convergence (optional, but good practice)
        if abs(theta_1_new - theta_1) < 1e-6 and abs(theta_2_new - theta_2) < 1e-6:
            break
            
        theta_1 = theta_1_new
        theta_2 = theta_2_new

    # 4. Calculate the final average number of layers
    avg_k = theta_1 + theta_2
    
    # Print the final equation with the calculated numbers
    print(f"The converged layer coverages are θ₁ ≈ {theta_1:.4f} and θ₂ ≈ {theta_2:.4f}.")
    print("The average number of adsorbed layers per site, <k>, is given by the sum:")
    print(f"<k> = θ₁ + θ₂ = {theta_1:.4f} + {theta_2:.4f} = {avg_k:.4f}")

    # Return the final numerical value as a string for the required output format
    return f"<<<{avg_k:.4f}>>>"

# Execute the function and print the result
final_answer = solve_adsorption()
print(final_answer)
