import math

def solve_adsorption():
    """
    Calculates the average number of adsorbed layers per site using a
    self-consistent mean-field model for a finite number of layers.
    """
    # Parameters from the problem (in units of k_B*T)
    beta_eps_1 = 0.1
    beta_mu = 0.15
    z_l = 4
    z_inter = 4

    # --- Assumptions to make the problem solvable ---
    # Assume the maximum number of layers 'k' is equal to the lateral coordination number z_l.
    k_max = z_l
    
    # Assume the interaction energy for subsequent layers (2 to k) is derived from eps_1 and z_inter.
    beta_eps_L = beta_eps_1 / z_inter
    
    # The lateral interaction energy depends on k_max.
    beta_eps_l = (0.02)**k_max

    # --- Self-consistent calculation for <k> ---
    # Initial guess for the average number of layers <k>
    avg_k = 1.0
    
    # Iterate to find the self-consistent solution for avg_k
    for _ in range(20): # 20 iterations are more than enough for convergence
        # Mean-field correction to the chemical potential (attractive interaction)
        # mu_eff = mu + U_mf, where U_mf = z_l * eps_l * <k>
        beta_mu_eff = beta_mu + z_l * beta_eps_l * avg_k
        
        # Terms in the partition function: exp(-beta(E_n - n*mu_eff))
        # Let T_n be the term for a stack of n particles
        # T_n = exp(beta*(eps_1 + (n-1)*eps_L + n*mu_eff)) for n>=1
        # Let x = exp(beta*(eps_1 + mu_eff)) and y = exp(beta*(eps_L + mu_eff))
        # T_n = x * y^(n-1)
        
        # Calculate x and y based on the current avg_k
        x = math.exp(beta_eps_1 + beta_mu_eff)
        y = math.exp(beta_eps_L + beta_mu_eff)

        # Calculate the single-site partition function xi = sum(T_n) for n=0 to k_max
        # T_0 = 1 (for an empty site)
        xi = 1.0
        # Calculate the numerator for <k>, num = sum(n * T_n)
        num = 0.0
        
        y_pow_n_minus_1 = 1.0
        for n in range(1, k_max + 1):
            term = x * y_pow_n_minus_1
            xi += term
            num += n * term
            y_pow_n_minus_1 *= y
            
        # Update the value of avg_k
        avg_k = num / xi

    # --- Output the final equation and result ---
    print("Based on the assumptions, the problem is solved by finding <k> in the equation:")
    print("<k> = (Σ n*T_n) / (Σ T_n), for n=0..k_max")
    print(f"where k_max = {k_max}\n")

    # Recalculate final terms for printing
    beta_mu_eff = beta_mu + z_l * beta_eps_l * avg_k
    x = math.exp(beta_eps_1 + beta_mu_eff)
    y = math.exp(beta_eps_L + beta_mu_eff)
    
    terms = [1.0]
    y_pow_n_minus_1 = 1.0
    for n in range(1, k_max + 1):
        terms.append(x * y_pow_n_minus_1)
        y_pow_n_minus_1 *= y

    numerator_str = " + ".join([f"{n}*{terms[n]:.4f}" for n in range(1, k_max + 1)])
    denominator_str = " + ".join([f"{t:.4f}" for t in terms])
    
    print("Final equation with calculated values:")
    print(f"<k> = ({numerator_str}) / ({denominator_str})")
    
    numerator_val = sum(n * terms[n] for n in range(1, k_max + 1))
    denominator_val = sum(terms)
    
    print(f"<k> = {numerator_val:.4f} / {denominator_val:.4f}")
    
    final_avg_k = numerator_val / denominator_val
    print(f"\nThe average number of adsorbed layers per site is <k> = {final_avg_k:.4f}")
    print(f"\n<<<${final_avg_k:.4f}>>>")

solve_adsorption()