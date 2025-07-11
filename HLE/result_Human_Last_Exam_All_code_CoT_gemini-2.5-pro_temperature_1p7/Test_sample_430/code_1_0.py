import math

def solve_adsorption():
    """
    Solves for the average number of adsorbed layers per site using a mean-field model.
    """
    # Parameters from the problem statement (in units of kB*T)
    beta_epsilon_1 = 0.1
    beta_mu = 0.15
    z_l = 4

    # Assumptions to make the problem solvable
    k_max_layers = 2 # Assumption for the max number of layers
    # Calculate beta*epsilon_l based on k=2
    beta_epsilon_l = (0.02)**k_max_layers
    # Assumption for vertical interaction energy (isotropic)
    beta_epsilon_v = beta_epsilon_l

    # --- Self-consistency equation: theta = G(theta) ---
    # For k=2, the equation for the average number of layers theta is:
    # theta = (Z1 + 2*Z2) / (1 + Z1 + Z2)
    # where Zm is the Boltzmann factor for a stack of m layers.
    # Zm = exp(beta*(m*mu - E_m(theta)))
    # E_m(theta) = epsilon_1 + (m-1)*epsilon_v + m*z_l*epsilon_l*theta

    # Let's define the terms in the exponents:
    # A(theta) = exp(beta*(mu - epsilon_1 - z_l*epsilon_l*theta))
    # x(theta) = exp(beta*(mu - epsilon_v - z_l*epsilon_l*theta))
    # Z1 = A(theta)
    # Z2 = A(theta) * x(theta)
    # So, theta = (A*(1 + 2x)) / (1 + A*(1+x))

    c1 = beta_mu - beta_epsilon_1
    c2 = z_l * beta_epsilon_l
    c3 = beta_mu - beta_epsilon_v

    # Iteratively solve for theta
    theta = 1.0  # Initial guess for the average number of layers
    for _ in range(20): # 20 iterations are more than enough for convergence
        A = math.exp(c1 - c2 * theta)
        x = math.exp(c3 - c2 * theta)
        
        numerator = A * (1 + 2 * x)
        denominator = 1 + A * (1 + x)
        
        theta_new = numerator / denominator
        
        # Check for convergence
        if abs(theta - theta_new) < 1e-9:
            theta = theta_new
            break
        theta = theta_new
        
    print("Based on the assumptions k=2 and epsilon_v = epsilon_l, we solve the self-consistency equation:")
    final_equation_str = (
        f"<L> = [ A * (1 + 2*x) ] / [ 1 + A * (1 + x) ]\n"
        f"where:\n"
        f"  A = exp(beta*(mu - e1) - beta*zl*el*<L>) = exp({c1:.4f} - {c2:.4f}*<L>)\n"
        f"  x = exp(beta*(mu - ev) - beta*zl*el*<L>) = exp({c3:.4f} - {c2:.4f}*<L>)\n"
    )
    print(final_equation_str)

    print(f"The calculated values for A and x at convergence are:")
    print(f"A = {A:.4f}")
    print(f"x = {x:.4f}\n")
    
    print("The final converged value for the average number of layers per site is:")
    print(f"<<<{theta:.4f}>>>")

solve_adsorption()