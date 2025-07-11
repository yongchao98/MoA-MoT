import numpy as np

def solve_adsorption_model():
    """
    Solves the self-consistent mean-field model for multi-layer adsorption.
    """
    # 1. Define constants from the problem statement (in units of k_B*T)
    beta_eps1 = 0.1
    beta_mu = 0.15
    z_l = 4.0
    z_inter = 4.0
    
    # Model assumptions
    k_max = 20  # Maximum number of layers to consider

    # 2. Iterative solver parameters
    max_iter_outer = 100
    max_iter_inner = 100
    tolerance = 1e-7
    mix_factor = 0.2  # Mixing factor for stability

    # 3. Outer self-consistency loop for <n>
    # Initial guess for the average number of layers, <n>
    n_avg = 1.0 
    
    final_thetas = np.zeros(k_max)

    for i in range(max_iter_outer):
        n_avg_old = n_avg

        # Calculate interaction energies based on current <n>
        if n_avg > 0:
            beta_eps_l = (0.02)**n_avg
        else:
            beta_eps_l = (0.02)**0 # Avoid math domain error for n_avg=0

        beta_eps_inter = beta_eps_l
        beta_eps2 = beta_eps_inter

        # 4. Inner self-consistency loop for {theta_i}
        # Initial guess for layer coverages
        thetas = np.full(k_max, 0.5)

        for j in range(max_iter_inner):
            thetas_old = np.copy(thetas)

            # Calculate effective energies H_i for each layer
            H_beta = np.zeros(k_max)
            H_beta[0] = -beta_eps1 - z_l * beta_eps_l * thetas[0]
            for l in range(1, k_max):
                H_beta[l] = -beta_eps2 - z_inter * beta_eps_inter - z_l * beta_eps_l * thetas[l]
            
            # Calculate total energy E_m for a site with m layers
            E_m_beta = np.cumsum(H_beta)

            # Calculate terms for the site partition function
            # exp(-beta*(E_m - m*mu)) = exp(m*beta*mu - beta*E_m)
            terms = np.zeros(k_max + 1)
            terms[0] = 1.0  # Weight for m=0 layers
            for m in range(1, k_max + 1):
                terms[m] = np.exp(m * beta_mu - E_m_beta[m-1])
            
            # Site partition function
            z_site = np.sum(terms)

            # Probability P(m) of having m layers
            P_m = terms / z_site
            
            # Update layer coverages theta_i = sum_{j=i to k_max} P(j)
            thetas_new = np.zeros(k_max)
            for l in range(k_max):
                thetas_new[l] = np.sum(P_m[l+1:])

            thetas = mix_factor * thetas_new + (1 - mix_factor) * thetas
            
            # Check for convergence of the inner loop
            if np.linalg.norm(thetas - thetas_old) < tolerance:
                break
        
        final_thetas = thetas
        n_avg_new = np.sum(thetas)
        
        # Update <n> with mixing
        n_avg = mix_factor * n_avg_new + (1 - mix_factor) * n_avg

        # Check for convergence of the outer loop
        if abs(n_avg - n_avg_old) < tolerance:
            break

    # 5. Output the results
    print("--- Model Parameters and Self-Consistent Solution ---")
    print("The average number of layers, <k>, is found by solving the following system of equations self-consistently.")
    print("\nUnderlying parameter values:")
    print(f"  beta * epsilon_1 = {beta_eps1}")
    print(f"  beta * mu = {beta_mu}")
    print(f"  z_l = {z_l}, z_inter = {z_inter}")

    print("\nSelf-consistent relationships (assumed):")
    print("  beta*epsilon_l = (0.02)^<k>")
    print("  beta*epsilon_2 = beta*epsilon_inter = beta*epsilon_l")
    
    print("\nConverged values:")
    final_beta_eps_l = (0.02)**n_avg
    print(f"  <k> = {n_avg:.8f}")
    print(f"  beta*epsilon_l = (0.02)^{n_avg:.4f} = {final_beta_eps_l:.8f}")

    H_beta_final = np.zeros(k_max)
    H_beta_final[0] = -beta_eps1 - z_l * final_beta_eps_l * final_thetas[0]
    for l in range(1, k_max):
        H_beta_final[l] = -final_beta_eps_l - z_inter * final_beta_eps_l - z_l * final_beta_eps_l * final_thetas[l]
    
    print("\nFinal equation is based on these converged effective energy terms (beta*H_i):")
    print(f"  beta*H_1 = -{beta_eps1} - {z_l}*{final_beta_eps_l:.6f}*{final_thetas[0]:.6f} = {H_beta_final[0]:.6f}")
    print(f"  beta*H_2 = -{final_beta_eps_l:.6f} - {z_inter}*{final_beta_eps_l:.6f} - {z_l}*{final_beta_eps_l:.6f}*{final_thetas[1]:.6f} = {H_beta_final[1]:.6f}")
    print("  ... and so on for higher layers.")
    
    print("\n--- Final Answer ---")
    print("The average number of adsorbed layers per site is:")
    print(n_avg)
    
    return n_avg

if __name__ == '__main__':
    average_layers = solve_adsorption_model()
    # The final answer format for the platform
    # print(f"<<<{average_layers:.8f}>>>")