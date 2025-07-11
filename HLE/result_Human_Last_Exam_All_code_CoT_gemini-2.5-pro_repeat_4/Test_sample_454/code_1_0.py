import math

def solve_lattice_gas_occupancy():
    """
    Calculates the average occupancy per site for a lattice gas system
    using the mean-field approximation and a fixed-point iteration method.
    """
    # 1. Define the parameters from the problem statement
    z_horiz = 4
    z_vert = 8
    mu_factor = 0.1  # mu = 0.1 * k_B * T
    
    # Total coordination number for a site in a 'bulk' layer
    z = z_horiz + z_vert
    
    # Dimensionless interaction energy: beta_epsilon = epsilon / (k_B * T)
    # epsilon = -k_B * T / (2 * pi)
    beta_epsilon = -1 / (2 * math.pi)
    
    # Dimensionless chemical potential: beta_mu = mu / (k_B * T)
    beta_mu = mu_factor

    print("Step 1: Theoretical Model")
    print("-------------------------")
    print("In the mean-field approximation, the average occupancy <n> is given by the self-consistency equation:")
    print("<n> = 1 / (exp(beta * (z * epsilon * <n> - mu)) + 1)")
    print("where beta = 1 / (k_B * T).\n")

    print("Step 2: Substitute Parameters")
    print("-----------------------------")
    print("We can rewrite the equation using dimensionless parameters:")
    print("<n> = 1 / (exp(z * (epsilon / (k_B * T)) * <n> - (mu / (k_B * T))) + 1)\n")
    print("The given parameters are:")
    print(f"Horizontal coordination number z_horizontial = {z_horiz}")
    print(f"Vertical coordination number z_vertical = {z_vert}")
    print(f"Total coordination number z = {z_horiz} + {z_vert} = {z}")
    print(f"Dimensionless interaction energy (epsilon / k_B*T) = -1 / (2*pi) = {beta_epsilon:.4f}")
    print(f"Dimensionless chemical potential (mu / k_B*T) = {beta_mu}\n")
    
    print("Step 3: Equation to Solve")
    print("--------------------------")
    print("Substituting these values, we need to solve the following equation for <n>:")
    print(f"<n> = 1 / (exp({z} * ({beta_epsilon:.4f}) * <n> - {beta_mu}) + 1)")
    print(f"<n> = 1 / (exp({z * beta_epsilon:.4f} * <n> - {beta_mu}) + 1)\n")

    # 2. Solve the equation using fixed-point iteration
    # Let g(n) be the right-hand side of the equation. We iterate n_new = g(n_old).
    
    # Initial guess for <n> (must be between 0 and 1)
    n = 0.5 
    
    # Number of iterations for convergence
    iterations = 20
    
    for _ in range(iterations):
        exponent = z * beta_epsilon * n - beta_mu
        n = 1 / (math.exp(exponent) + 1)
        
    # 3. Print the final result
    print("Step 4: Numerical Solution")
    print("--------------------------")
    print("Solving the equation numerically using fixed-point iteration yields:")
    print(f"The average occupancy per site <n> is: {n:.3f}")
    
    # Final answer in the required format
    print(f"\n<<<{n:.3f}>>>")

if __name__ == "__main__":
    solve_lattice_gas_occupancy()