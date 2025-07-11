import numpy as np

def solve_adsorption():
    """
    Solves for the average number of adsorbed layers per site using a
    self-consistent mean-field model based on the provided parameters.
    """
    # Parameters from the problem statement
    beta_eps1 = 0.1  # beta * epsilon_1
    beta_mu = 0.15   # beta * mu
    zl = 4           # lateral coordination number
    z_inter = 4      # vertical coordination number

    # --- Assumptions to resolve ambiguities in the problem statement ---
    # 1. Assume the maximum number of layers, k, is 4 (same as coordination numbers).
    k_max = 4
    # 2. Assume eps_inter = eps_l, and eps_2 = z_inter * eps_inter.
    # This makes all interaction energies depend on k_max.
    # 3. The mean-field lateral interaction energy for a j-stack is j * zl * eps_l * <k>.

    # Calculate derived energy parameters (in units of k_B*T)
    beta_eps_l = (0.02)**k_max
    beta_eps2 = z_inter * beta_eps_l
    
    # Pre-calculate terms for the self-consistent equation
    # These are the layer activities without the mean-field term
    q0 = np.zeros(k_max + 1)
    q0[0] = 1.0 # for the empty state
    # First layer activity (without MF)
    q0[1] = np.exp(-(beta_eps1 - beta_mu))
    # Subsequent layer activities (without MF)
    for j in range(2, k_max + 1):
        # E_j = eps1 + (j-1)*eps2
        # E_j - j*mu = (eps1 - mu) + (j-1)*(eps2 - mu) - (j-1)*mu = WRONG
        # E_j - j*mu = (eps1 - mu) + (j-1)*eps2 - (j-1)*mu
        # beta*(E_j - j*mu) = beta*(eps1-mu) + (j-1)*beta*eps2 - (j-1)*beta*mu
        # = beta*(eps1-mu) + (j-1)*beta*(eps2-mu)
        beta_Ej_minus_jmu = (beta_eps1 - beta_mu) + (j - 1) * (beta_eps2 - beta_mu)
        q0[j] = np.exp(-beta_Ej_minus_jmu)
        
    # Self-consistent calculation for n = <k> (average number of layers)
    # Using fixed-point iteration: n_new = F(n_old)
    
    n = 1.0  # Initial guess for the average number of layers
    for _ in range(100): # Iterate to find a stable solution
        # Mean-field term exponent factor
        mf_factor = zl * beta_eps_l
        
        # Calculate partition function terms q_j with mean-field correction
        q = np.zeros(k_max + 1)
        q[0] = 1.0
        for j in range(1, k_max + 1):
            # E_eff = E_vert + j*zl*eps_l*n - j*mu
            # beta*E_eff = beta*(E_vert - j*mu) + j*zl*beta_eps_l*n
            # exp(-beta*E_eff) = q0[j] * exp(-j*mf_factor*n)
            q[j] = q0[j] * np.exp(-j * mf_factor * n)

        # Calculate the numerator and denominator for <k>
        numerator = np.sum([j * q[j] for j in range(1, k_max + 1)])
        denominator = np.sum(q)
        
        n_new = numerator / denominator
        
        if np.abs(n_new - n) < 1e-9: # Convergence check
            break
        n = n_new

    # --- Output the results ---
    # The final average number of layers
    print(f"The average number of adsorbed layers per site is: {n:.6f}")
    
    # The final self-consistent equation with numerical values
    print("\nThis value is the solution 'n' to the self-consistent equation: n = F(n)")
    print("F(n) = (Σ_{j=1 to k_max} j * q_j(n)) / (Σ_{j=0 to k_max} q_j(n))\n")
    print(f"Based on the assumptions (k_max={k_max}), the equation with the final value of n={n:.6f} is:")
    
    mf_factor_final = zl * beta_eps_l
    q_final = np.zeros(k_max + 1)
    q_final[0] = 1.0
    for j in range(1, k_max + 1):
        q_final[j] = q0[j] * np.exp(-j * mf_factor_final * n)
        
    numerator_str = " + ".join([f"{j}*exp(-({beta_eps1 - beta_mu:.4f} + {(j-1)*(beta_eps2 - beta_mu):.4f} + {j*mf_factor_final*n:.6f}))" for j in range(1, k_max + 1)])
    denominator_str = "1.0 + " + " + ".join([f"exp(-({beta_eps1 - beta_mu:.4f} + {(j-1)*(beta_eps2 - beta_mu):.4f} + {j*mf_factor_final*n:.6f}))" for j in range(1, k_max + 1)])
    
    numerator_val = np.sum([j * q_final[j] for j in range(1, k_max + 1)])
    denominator_val = np.sum(q_final)

    print(f"n = ({numerator_val:.6f}) / ({denominator_val:.6f})")

    # The problem asks for each number in the final equation. 
    # Printing the components of the sum more explicitly.
    print("\nWhere the numerator is the sum of j*q_j(n):")
    num_terms = []
    for j in range(1, k_max + 1):
        term_val = j * q_final[j]
        num_terms.append(f"{term_val:.6f}")
    print(" + ".join(num_terms))

    print("\nAnd the denominator is the sum of q_j(n):")
    den_terms = ["1.0"]
    for j in range(1, k_max + 1):
        term_val = q_final[j]
        den_terms.append(f"{term_val:.6f}")
    print(" + ".join(den_terms))
    print(f"\n<<< {n:.6f} >>>")

solve_adsorption()