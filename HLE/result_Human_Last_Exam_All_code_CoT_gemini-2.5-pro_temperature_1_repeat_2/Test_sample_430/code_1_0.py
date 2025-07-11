import numpy as np

def solve_adsorption_model():
    """
    Solves the multilayer adsorption problem using a self-consistent mean-field model.
    """
    # 1. Define parameters based on the problem statement and interpretation
    # Energies are in units of k_B*T, so we can set k_B*T = 1
    mu = 0.15          # Chemical potential
    eps1_vert = 0.1    # Vertical binding energy for the first layer
    z_l = 4            # Lateral coordination number
    
    # Interpreted parameters based on the plan
    k_max = 20         # Maximum number of layers to consider
    eps_l = 0.02       # Assumed lateral interaction energy
    z_inter = 4        # Vertical coordination factor
    
    # Derived vertical energy for layers s >= 2
    # eps_s_vert = z_inter * eps_l
    eps_s_vert = z_inter * eps_l 
    
    # Pre-calculate the mean-field lateral interaction term factor
    zl_epsl = z_l * eps_l

    # 2. Iteratively solve the self-consistent equations
    max_iter = 1000
    tolerance = 1e-9
    
    # Initial guess for layer coverages theta_s
    # theta[s] corresponds to theta_{s+1}
    theta = np.zeros(k_max)

    for i in range(max_iter):
        # Calculate the ratios x_s / x_{s-1}
        # R[s] corresponds to ratio for x_{s+1}
        R = np.zeros(k_max)
        
        # For layer s=1
        exponent1 = mu + eps1_vert + zl_epsl * theta[0]
        R[0] = np.exp(exponent1)
        
        # For layers s=2 to k_max
        for s in range(1, k_max):
            exponent_s = mu + eps_s_vert + zl_epsl * theta[s]
            R[s] = np.exp(exponent_s)

        # Calculate unnormalized probabilities x_s
        x_unnormalized = np.zeros(k_max + 1)
        x_unnormalized[0] = 1.0
        for s in range(1, k_max + 1):
            x_unnormalized[s] = x_unnormalized[s-1] * R[s-1]
            
        # Normalize probabilities
        total_x = np.sum(x_unnormalized)
        x = x_unnormalized / total_x
        
        # Calculate new theta values from the new x values
        theta_new = np.zeros(k_max)
        for s in range(k_max):
            theta_new[s] = np.sum(x[s+1:])
            
        # Check for convergence
        if np.sum(np.abs(theta_new - theta)) < tolerance:
            break
            
        theta = theta_new

    # 3. Calculate the final average number of layers, <k>
    # <k> = sum(s * x_s) for s from 1 to k_max
    avg_k = 0.0
    calculation_str_parts = []
    for s in range(1, k_max + 1):
        avg_k += s * x[s]
        if x[s] > 1e-6: # Only print significant contributions
             calculation_str_parts.append(f"({s} * {x[s]:.4f})")

    # 4. Print the final result and the calculation
    print("Based on the self-consistent mean-field model, the converged probabilities for a site having 's' layers (x_s) are:")
    for s in range(k_max + 1):
        if x[s] > 1e-6:
            print(f"x_{s:<2} = {x[s]:.6f}")

    print("\nThe average number of adsorbed layers per site, <k>, is calculated as:")
    calculation_str = " <k> = " + " + ".join(calculation_str_parts)
    print(calculation_str)
    print(f"\n<k> = {avg_k}")
    
    # Final answer format
    print(f"\n<<<{avg_k}>>>")


solve_adsorption_model()