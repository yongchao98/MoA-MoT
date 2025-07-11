import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site for a lattice gas
    using the mean-field approximation.
    """
    # Given parameters
    # beta_epsilon = (epsilon) / (k_B * T)
    beta_epsilon = -1 / (2 * math.pi)
    
    # beta_mu = mu / (k_B * T)
    beta_mu = 0.1
    
    # Coordination numbers
    z_horizontal = 4
    z_vertical = 8
    z = z_horizontal + z_vertical
    
    # The self-consistency equation is: <n> = 1 / (1 + exp(A * <n> + B))
    # where A = z * beta_epsilon and B = -beta_mu
    A = z * beta_epsilon
    B = -beta_mu

    print("The self-consistency equation is of the form: <n> = 1 / (1 + exp(A * <n> + B))")
    print(f"The numbers in the final equation are:")
    print(f"A = z * (epsilon / k_B*T) = {z} * ({beta_epsilon:.4f}) = {A:.4f}")
    print(f"B = -mu / (k_B*T) = {-beta_mu:.4f}")
    print("-" * 20)

    # Solve for <n> using fixed-point iteration
    n_avg = 0.5  # Initial guess
    
    # Iterate to find the stable solution
    for _ in range(100):
        exponent = A * n_avg + B
        n_avg_new = 1 / (1 + math.exp(exponent))
        
        # Check for convergence
        if abs(n_avg_new - n_avg) < 1e-6:
            break
        n_avg = n_avg_new

    print(f"The calculated average occupancy per site <n> is: {n_avg:.3f}")

solve_average_occupancy()
<<<0.848>>>