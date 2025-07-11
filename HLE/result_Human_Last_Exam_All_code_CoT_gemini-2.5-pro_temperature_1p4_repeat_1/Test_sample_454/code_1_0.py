import numpy as np

def solve_lattice_gas_occupancy():
    """
    Calculates the average site occupancy for a lattice gas model using the
    mean-field approximation.
    """
    # Step 1: Define the parameters based on the problem statement.
    # We work in dimensionless units where beta = 1/(k_B*T) is absorbed into the energy terms.
    # beta * epsilon = (1/(k_B*T)) * (-k_B*T / (2*pi))
    beta_epsilon = -1.0 / (2.0 * np.pi)
    
    # beta * mu = (1/(k_B*T)) * (0.1 * k_B*T)
    beta_mu = 0.1
    
    z_horiz = 4
    z_vert = 8
    # The total coordination number z for a site not in the monolayer.
    z = z_horiz + z_vert
    
    # Step 2: Print the theoretical context and the equation being solved.
    print("We are solving the self-consistency equation for the average occupancy <n> in the mean-field approximation:")
    print("      <n> = 1 / (exp(beta*(epsilon*z*<n> - mu)) + 1)")
    print("\nWith the given parameters:")
    print(f"      beta*epsilon = -1/(2*pi) = {beta_epsilon:.4f}")
    print(f"      beta*mu = 0.1")
    print(f"      z = z_horizontal + z_vertical = {z_horiz} + {z_vert} = {z}")
    
    # Combined parameter for the exponent
    beta_epsilon_z = beta_epsilon * z
    
    print("\nThe specific equation to solve is:")
    # Print the equation with each number explicitly shown
    print(f"      <n> = 1 / (exp( ({beta_epsilon:.4f}) * ({z}) * <n> - ({beta_mu:.1f}) ) + 1)")
    print(f"      <n> = 1 / (exp( {beta_epsilon_z:.4f} * <n> - {beta_mu:.1f} ) + 1)")
    
    # Step 3: Solve the equation numerically using an iterative method.
    # Start with an initial guess for <n> (must be between 0 and 1)
    n = 0.5 
    # Iterate to find a stable solution
    for _ in range(100):
        exponent = beta_epsilon_z * n - beta_mu
        n_new = 1.0 / (np.exp(exponent) + 1.0)
        # Check for convergence
        if abs(n_new - n) < 1e-7:
            n = n_new
            break
        n = n_new
        
    # Step 4: Print the final result.
    print("\nAfter solving numerically, the average occupancy per site is:")
    print(f"      <n> = {n:.3f}")
    
    # Final answer in the required format
    print(f"<<<{n:.3f}>>>")

solve_lattice_gas_occupancy()