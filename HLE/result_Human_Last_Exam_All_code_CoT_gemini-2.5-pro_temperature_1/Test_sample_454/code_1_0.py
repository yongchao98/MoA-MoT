import numpy as np

def solve_average_occupancy():
    """
    Calculates the average occupancy per site <n> for a lattice gas
    using the mean-field approximation and fixed-point iteration.
    """
    # Given parameters
    z_horizontal = 4
    z_vertical = 8
    
    # Total coordination number z
    z = z_horizontal + z_vertical
    
    # Dimensionless parameters based on the problem statement
    # mu = 0.1 * k_B * T  => beta * mu = 0.1
    # epsilon = -k_B * T / (2 * pi) => beta * epsilon = -1 / (2 * pi)
    beta_mu = 0.1
    beta_epsilon = -1 / (2 * np.pi)

    # The self-consistency equation is:
    # <n> = 1 / (1 + exp(-(beta*mu + z * beta*epsilon * <n>)))
    
    print("The self-consistency equation for the average occupancy <n> is:")
    print(f"<n> = 1 / (1 + exp(-({beta_mu} + {z} * ({beta_epsilon:.4f}) * <n>)))")
    print(f"<n> = 1 / (1 + exp(-({beta_mu} - {z/(2*np.pi):.4f} * <n>)))")
    print("\nSolving this equation numerically...\n")

    # --- Numerical Solution using Fixed-Point Iteration ---
    n_avg = 0.5  # Initial guess for <n>
    tolerance = 1e-7
    
    for i in range(100): # Iterate to find the solution
        exponent = -(beta_mu + z * beta_epsilon * n_avg)
        n_avg_next = 1.0 / (1.0 + np.exp(exponent))
        
        # Check for convergence
        if abs(n_avg_next - n_avg) < tolerance:
            break
        
        n_avg = n_avg_next

    # Final result
    print(f"The calculated average occupancy per site <n> is: {n_avg:.3f}")
    return n_avg

# Execute the function
final_n_avg = solve_average_occupancy()
# The final answer is wrapped as requested
# print(f"<<<{final_n_avg:.3f}>>>")