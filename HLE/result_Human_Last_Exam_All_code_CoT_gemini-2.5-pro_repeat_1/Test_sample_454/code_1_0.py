import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site <n> for a lattice gas
    using the mean-field approximation.
    """
    # Given parameters
    mu_over_kBT = 0.1
    epsilon_over_kBT = -1 / (2 * math.pi)
    z_horizontal = 4
    z_vertical = 8
    
    # In the mean-field model assuming all sites are equivalent, the total
    # coordination number is the sum of horizontal and vertical neighbors.
    z = z_horizontal + z_vertical
    
    # The self-consistency equation for the average occupancy <n> (theta) is:
    # theta = 1 / (1 + exp(z * (epsilon/kBT) * theta - (mu/kBT)))
    
    # --- Numerical Solution using Fixed-Point Iteration ---
    # Initial guess for the average occupancy <n>
    n_avg = 0.5 
    
    # Iterate to find the solution. 100 iterations are sufficient for convergence.
    for _ in range(100):
        exponent = z * epsilon_over_kBT * n_avg - mu_over_kBT
        n_avg = 1 / (1 + math.exp(exponent))

    # --- Output the Results ---
    print("This script solves for the average occupancy <n> in a lattice gas model.")
    print("The self-consistency equation derived from mean-field theory is:")
    print("<n> = 1 / (1 + exp( (z * epsilon / (kB * T)) * <n> - (mu / (kB * T)) ))\n")
    
    print("Substituting the given parameter values:")
    print(f"z = {z_horizontal} + {z_vertical} = {z}")
    print(f"mu / (kB * T) = {mu_over_kBT}")
    print(f"epsilon / (kB * T) = -1 / (2*pi) approx {epsilon_over_kBT:.4f}\n")
    
    z_eps_term_str = f"-{z}/(2*pi)"
    if z % 2 == 0:
      z_eps_term_str = f"-{z//2}/pi"

    print("The final equation to be solved for <n> is:")
    print(f"<n> = 1 / (1 + exp( ({z_eps_term_str}) * <n> - {mu_over_kBT} ))\n")
    
    # Print the final calculated value rounded to three decimal places
    print(f"The calculated average occupancy per site <n> is: {n_avg:.3f}")

solve_average_occupancy()
<<<0.848>>>