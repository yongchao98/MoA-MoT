import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site for a lattice gas model
    using the mean-field approximation.
    """
    # Given parameters
    z_horizontial = 4
    z_vertical = 8
    mu_over_kBT = 0.1
    
    # The temperature T=300K is not explicitly needed as epsilon and mu are given
    # in terms of k_B*T.
    # beta * epsilon = (1/(k_B*T)) * (-k_B*T / (2*pi))
    epsilon_over_kBT = -1 / (2 * math.pi)
    
    # Total coordination number
    z = z_horizontial + z_vertical
    
    # --- Numerical Iterative Solver ---
    # Initial guess for the average occupancy <n>
    n_avg = 0.5
    
    # Iteration parameters
    tolerance = 1e-7
    max_iterations = 100
    
    print("Solving the self-consistency equation for the average occupancy <n>:")
    print("<n> = 1 / (exp(z * (ε/kBT) * <n> - (μ/kBT)) + 1)\n")
    
    print("With the following parameters:")
    print(f"  Total coordination number, z = {z_horizontial} + {z_vertical} = {z}")
    print(f"  Interaction energy term, ε/kBT = -1/(2π) ≈ {epsilon_over_kBT:.4f}")
    print(f"  Chemical potential term, μ/kBT = {mu_over_kBT:.4f}\n")
    
    print("The final equation to solve is:")
    final_coeff = z * epsilon_over_kBT
    print(f"<n> = 1 / (exp({z} * {epsilon_over_kBT:.4f} * <n> - {mu_over_kBT:.4f}) + 1)")
    print(f"<n> = 1 / (exp({final_coeff:.4f} * <n> - {mu_over_kBT:.4f}) + 1)\n")
    
    for i in range(max_iterations):
        # Update <n> using the self-consistency equation
        exponent = z * epsilon_over_kBT * n_avg - mu_over_kBT
        n_avg_new = 1 / (math.exp(exponent) + 1)
        
        # Check for convergence
        if abs(n_avg_new - n_avg) < tolerance:
            n_avg = n_avg_new
            break
        
        n_avg = n_avg_new
        
    print(f"After {i+1} iterations, the solution converges.")
    print(f"The calculated average occupancy per site <n> is: {n_avg:.3f}")

    return n_avg

# Run the solver
final_occupancy = solve_average_occupancy()