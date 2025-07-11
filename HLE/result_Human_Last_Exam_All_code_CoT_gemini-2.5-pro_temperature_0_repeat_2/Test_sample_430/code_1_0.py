import math

def solve_average_layers():
    """
    Calculates the average number of adsorbed layers per site, <k>,
    by solving a self-consistent mean-field equation.
    """
    # Parameters from the problem, normalized by k_B*T
    mu_prime = 0.15
    eps1_prime = 0.1
    z_l = 4
    z_inter = 4

    # Initial guess for <k>, the average number of layers
    k_avg = 1.0
    
    # Use fixed-point iteration to find the self-consistent value of <k>
    for i in range(100):  # 100 iterations are sufficient for convergence
        k_old = k_avg
        
        # 1. Calculate the lateral interaction energy, which depends on <k>
        # ε'_ℓ = (0.02)^<k>
        eps_l_prime = (0.02)**k_avg
        
        # 2. Assume vertical energy for subsequent layers is ε'_s = z_inter * ε'_ℓ
        eps_s_prime = z_inter * eps_l_prime
        
        # 3. Apply the simplified mean-field approximation for the interaction term.
        # The interaction energy is approximated using the average <k>
        mf_term = z_l * eps_l_prime * k_avg
        
        # 4. Calculate effective energies for the first and subsequent layers
        eff_eps1_prime = eps1_prime + mf_term
        eff_eps_s_prime = eps_s_prime + mf_term
        
        # 5. Calculate the adsorption parameters x1 and xs, analogous to the BET model
        # x1 = exp(μ' - ε'_eff,1)
        x1 = math.exp(mu_prime - eff_eps1_prime)
        
        # xs = exp(μ' - ε'_eff,s)
        xs = math.exp(mu_prime - eff_eps_s_prime)
        
        # 6. Update <k> using the derived formula <k> = G(<k>)
        # This formula gives the total number of layers as a sum of a geometric series of layer occupancies
        if xs >= 1:
            # This case corresponds to bulk condensation, <k> would be infinite.
            # Based on the parameters, this should not happen.
            k_avg = float('inf')
            break
        
        k_avg = (x1 / (1 + x1)) * (1 + xs)
        
        # 7. Check for convergence
        if abs(k_avg - k_old) < 1e-9:
            break
            
    # Output the results as requested
    print("The self-consistent equation to be solved is of the form: <k> = G(<k>)")
    print("where G(<k>) = (x1 / (1 + x1)) * (1 + xs), and x1, xs depend on <k>.")
    print("\nAfter convergence, the final equation with the calculated values is:")
    print(f"{k_avg:.5f} = ({x1:.5f} / (1 + {x1:.5f})) * (1 + {xs:.5f})")
    
    print("\nThe average number of adsorbed layers per site, <k>, is:")
    print(k_avg)
    
    return k_avg

# Run the calculation and store the final answer
final_k_value = solve_average_layers()