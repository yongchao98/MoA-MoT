import math

def solve_adsorption_model():
    """
    Solves for the average number of adsorbed layers per site using a
    mean-field model for a bilayer system (k=2).
    """
    # Parameters from the problem statement (in units of k_B*T)
    beta_e1 = 0.1
    beta_mu = 0.15
    z_l = 4
    z_inter = 4
    k_max = 2  # Assumption: maximum number of layers is 2

    # Calculate derived parameters based on assumptions
    # beta_e_ell = beta * epsilon_ell
    beta_e_ell = (0.02)**k_max
    # beta_e2 = beta * epsilon_2, assuming e2 = z_inter * e_ell
    beta_e2 = z_inter * beta_e_ell

    # Self-consistent calculation for <L> (average layers)
    # Initial guess for L
    L = k_max / 2.0
    
    # Iterate to find the self-consistent solution for L
    for _ in range(100): # 100 iterations is more than enough for convergence
        theta = L / k_max
        
        # Grand potential for a site with 1 layer (in units of k_B*T)
        # beta_Omega_1 = beta*(e1 - mu) + beta*z_l*e_ell*theta
        beta_Omega_1 = (beta_e1 - beta_mu) + z_l * beta_e_ell * theta

        # Grand potential for a site with 2 layers (in units of k_B*T)
        # beta_Omega_2 = beta*(e1 + e2 - 2*mu) + beta*2*z_l*e_ell*theta
        beta_Omega_2 = (beta_e1 + beta_e2 - 2 * beta_mu) + 2 * z_l * beta_e_ell * theta

        # Probabilities (unnormalized)
        prob1 = math.exp(-beta_Omega_1)
        prob2 = math.exp(-beta_Omega_2)
        
        # Partition function for a single site
        Xi = 1 + prob1 + prob2
        
        # New value for L
        L_new = (1 * prob1 + 2 * prob2) / Xi
        
        # Check for convergence
        if abs(L - L_new) < 1e-9:
            L = L_new
            break
        L = L_new

    # Final calculation with converged L to display the equation
    theta_final = L / k_max
    beta_Omega_1_final = (beta_e1 - beta_mu) + z_l * beta_e_ell * theta_final
    beta_Omega_2_final = (beta_e1 + beta_e2 - 2 * beta_mu) + 2 * z_l * beta_e_ell * theta_final
    prob1_final = math.exp(-beta_Omega_1_final)
    prob2_final = math.exp(-beta_Omega_2_final)
    Xi_final = 1 + prob1_final + prob2_final
    
    print("Based on the mean-field model for a bilayer system (k=2):")
    print("\nThe self-consistent equation for the average number of layers <L> is:")
    print(f"<L> = (1 * exp(-beta*\u03A9_1) + 2 * exp(-beta*\u03A9_2)) / (1 + exp(-beta*\u03A9_1) + exp(-beta*\u03A9_2))")
    print("\nPlugging in the converged value <L> = {:.4f}:".format(L))
    print("{:.4f} = (1 * exp(-({:.4f})) + 2 * exp(-({:.4f}))) / (1 + exp(-({:.4f})) + exp(-({:.4f})))".format(
        L, beta_Omega_1_final, beta_Omega_2_final, beta_Omega_1_final, beta_Omega_2_final))
    print("{:.4f} = (1 * {:.4f} + 2 * {:.4f}) / (1 + {:.4f} + {:.4f})".format(
        L, prob1_final, prob2_final, prob1_final, prob2_final))
    print("{:.4f} = {:.4f} / {:.4f}".format(
        L, (prob1_final + 2 * prob2_final), Xi_final))

    print(f"\nThe average number of adsorbed layers per site is {L:.4f}")

    # Return the final answer in the required format
    return L

if __name__ == '__main__':
    average_layers = solve_adsorption_model()
    # The final answer format for the platform
    # print(f"\n<<<{average_layers:.4f}>>>")
    # For this specific problem, the final output needs to be printed directly. Let's find the result and format it.
    print(f"\n<<<{average_layers:.4f}>>>")
