import math

def solve_adsorption_occupancy():
    """
    Solves for the average occupancy in a lattice gas model using mean-field theory.
    """
    # 1. Define the given parameters in dimensionless form (divided by kB*T)
    mu_over_kBT = 0.1
    epsilon_over_kBT = -1.0 / (2.0 * math.pi)
    z_horizontial = 4
    z_vertical = 8
    
    # Total coordination number z is the sum of its components
    z = z_horizontial + z_vertical

    # 2. Set up the self-consistency equation for the average occupancy <n>, denoted as n.
    # The equation is: n = 1 / (1 + exp(-(mu' - z * epsilon' * n)))
    # where mu' = mu / (kB*T) and epsilon' = epsilon / (kB*T).
    
    # We can write this as: n = 1 / (1 + exp(C1 + C2 * n))
    C1 = -mu_over_kBT
    C2 = z * epsilon_over_kBT

    print("Based on the mean-field approximation, the self-consistency equation is:")
    print("<n> = 1 / (1 + exp(-(mu/(kB*T) - z*epsilon/(kB*T) * <n>)))")
    print("\nWith the given values, the equation becomes:")
    print(f"<n> = 1 / (1 + exp(-({mu_over_kBT} - {z} * ({epsilon_over_kBT:.5f}) * <n>)))")
    print(f"<n> = 1 / (1 + exp({C1:.5f} + {C2:.5f} * <n>))")
    
    # 3. Solve the equation numerically using fixed-point iteration
    n = 0.5  # Initial guess for the average occupancy
    
    # Iterate to find the solution until it converges
    for _ in range(50):
        n_new = 1.0 / (1.0 + math.exp(C1 + C2 * n))
        if abs(n_new - n) < 1e-7: # Convergence check
            n = n_new
            break
        n = n_new

    # 4. Round the final answer to three decimal places and print it.
    final_n = round(n, 3)

    print("\nSolving for <n> numerically gives:")
    print(f"The average occupancy per site <n> = {final_n}")


# Execute the solver function
solve_adsorption_occupancy()
<<<0.848>>>