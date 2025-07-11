import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site using the mean-field approximation.
    
    The self-consistency equation is:
    <n> = 1 / (1 + exp(-beta_mu + z * beta_epsilon * <n>))
    
    We solve this equation for <n> using an iterative method.
    """
    
    # Given parameters
    mu_over_kT = 0.1  # This is beta * mu
    epsilon_over_kT = -1 / (2 * math.pi) # This is beta * epsilon
    z_h = 4
    z_v = 8
    
    # Total coordination number in the mean-field approximation
    z = z_h + z_v
    
    # --- Numerical Solver for <n> ---
    # Initial guess for the average occupancy <n>
    n_avg = 0.5 
    
    # Iterate to find the self-consistent solution
    # 100 iterations are more than enough for convergence here
    for _ in range(100):
        exponent = -mu_over_kT + z * epsilon_over_kT * n_avg
        n_avg = 1 / (1 + math.exp(exponent))
        
    # --- Output the results ---
    
    # Final equation with numerical values plugged in
    equation_str = (f"<n> = 1 / (1 + exp(-{mu_over_kT} + "
                    f"({z_h} + {z_v}) * ({epsilon_over_kT:.4f}) * <n>))")

    print("The self-consistency equation to be solved is:")
    print(equation_str)
    print("\nThis simplifies to:")
    print(f"<n> = 1 / (1 + exp(-{mu_over_kT} - {abs(z * epsilon_over_kT):.4f} * <n>))")
    
    # Print the final calculated value
    print(f"\nThe calculated average occupancy per site <n> is:")
    print(f"{n_avg:.3f}")

# Execute the function to find and print the solution
solve_average_occupancy()
<<<0.848>>>