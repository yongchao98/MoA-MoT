import math

def solve_adsorption_occupancy():
    """
    Solves for the average occupancy per site in a lattice gas model
    using the mean-field approximation.
    """
    # --- Step 1: Define Parameters ---
    # The term k_B*T cancels out when we use dimensionless quantities beta*epsilon and beta*mu.
    # beta = 1 / (k_B * T)
    
    # Interaction energy term: epsilon = -k_B*T / (2*pi)
    beta_epsilon = -1 / (2 * math.pi)
    
    # Chemical potential term: mu = 0.1 * k_B*T
    beta_mu = 0.1
    
    # Coordination numbers
    z_horizontal = 4
    z_vertical = 8
    
    # Total coordination number for a bulk-like particle
    z = z_horizontal + z_vertical

    # --- Step 2: Set up and solve the self-consistent equation ---
    # The equation is: n = 1 / (exp(z * beta_epsilon * n - beta_mu) + 1)
    # We solve this using fixed-point iteration: n_new = f(n_old)

    # Initial guess for the average occupancy <n>
    n = 0.5 
    
    # Iteration parameters
    tolerance = 1e-7
    max_iterations = 100
    
    for i in range(max_iterations):
        n_old = n
        # This is the function f(n) we are iterating
        exponent = z * beta_epsilon * n_old - beta_mu
        n = 1 / (math.exp(exponent) + 1)
        
        # Check for convergence
        if abs(n - n_old) < tolerance:
            break
    
    # --- Step 3: Print the derivation and results ---
    print("The self-consistent mean-field equation for the average occupancy <n> is:")
    print("n = 1 / (exp(z * (epsilon / k_B*T) * n - (mu / k_B*T)) + 1)")
    print("\nSubstituting the given parameter values:")
    
    z_beta_epsilon_val = z * beta_epsilon
    
    print(f"z = z_horizontal + z_vertical = {z_horizontal} + {z_vertical} = {z}")
    print(f"epsilon / (k_B*T) = -1 / (2*pi) ≈ {beta_epsilon:.5f}")
    print(f"mu / (k_B*T) = {beta_mu}")
    
    print("\nThe equation becomes:")
    print(f"n = 1 / (exp({z} * ({beta_epsilon:.5f}) * n - {beta_mu}) + 1)")
    print(f"n = 1 / (exp({z_beta_epsilon_val:.5f} * n - {beta_mu}) + 1)")
    
    print("\nSolving this equation numerically yields:")
    print(f"Average occupancy per site <n> = {n:.3f}")

# Execute the function
solve_adsorption_occupancy()
<<<0.848>>>