import math

def solve_lattice_gas():
    """
    Calculates the average occupancy per site for a lattice gas model
    using the mean-field approximation.
    """
    # 1. Define the parameters from the problem statement
    z_horizontal = 4
    z_vertical = 8
    
    # We use a total coordination number z for a representative "bulk" particle
    z = z_horizontal + z_vertical 
    
    # The parameters are given in terms of k_B*T, so we can use beta*E directly
    # beta = 1 / (k_B * T)
    # beta * epsilon = (1/(k_B*T)) * (-k_B*T / (2*pi)) = -1/(2*pi)
    beta_epsilon = -1.0 / (2.0 * math.pi)
    
    # beta * mu = (1/(k_B*T)) * (0.1 * k_B*T) = 0.1
    beta_mu = 0.1
    
    # 2. Print the self-consistency equation being solved
    print("The self-consistency equation for the average occupancy <n> is:")
    print("<n> = 1 / (1 + exp(beta * (z * epsilon * <n> - mu)))")
    print("\nPlugging in the given values:")
    print(f"<n> = 1 / (1 + exp({z} * ({beta_epsilon:.4f}) * <n> - {beta_mu}))")
    
    z_beta_epsilon = z * beta_epsilon
    print(f"<n> = 1 / (1 + exp({z_beta_epsilon:.4f} * <n> - {beta_mu}))")

    # 3. Solve the transcendental equation using fixed-point iteration
    n = 0.5  # Initial guess for <n>
    tolerance = 1e-7
    max_iterations = 100
    
    for i in range(max_iterations):
        exponent = z_beta_epsilon * n - beta_mu
        n_new = 1.0 / (1.0 + math.exp(exponent))
        
        if abs(n_new - n) < tolerance:
            n = n_new
            break
        n = n_new
    else:
        # This part runs if the loop finishes without converging
        print("\nWarning: The calculation did not converge within the maximum iterations.")

    # 4. Print the final result
    print(f"\nSolved average occupancy per site <n> = {n:.3f}")

solve_lattice_gas()