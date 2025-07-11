import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site <n> for a lattice gas
    using the mean-field approximation and fixed-point iteration.
    """
    # Given parameters
    z_horizontal = 4
    z_vertical = 8
    
    # Total coordination number z is the sum of horizontal and vertical neighbors
    z = z_horizontal + z_vertical

    # The self-consistent equation for the average occupancy <n> is:
    # <n> = 1 / (1 + exp(-beta * (mu - z * epsilon * <n>)))
    #
    # We are given:
    # mu = 0.1 * k_B * T  => beta * mu = 0.1
    # epsilon = -k_B * T / (2 * pi) => beta * epsilon = -1 / (2 * pi)
    #
    # Substituting these gives:
    # <n> = 1 / (1 + exp(-(0.1 - (-1 / (2 * pi)) * z * <n>)))
    # <n> = 1 / (1 + exp(-(0.1 + (z / (2 * pi)) * <n>)))
    
    # We solve this equation for <n> using fixed-point iteration.
    # Let's define the function f(<n>) = RHS. We iterate <n>_new = f(<n>_old).
    
    n_avg = 0.5  # Initial guess for <n>
    tolerance = 1e-9
    max_iterations = 1000

    beta_mu = 0.1
    beta_epsilon_z_term = z / (2 * math.pi)

    for i in range(max_iterations):
        exponent = -(beta_mu + beta_epsilon_z_term * n_avg)
        n_avg_new = 1 / (1 + math.exp(exponent))
        
        if abs(n_avg_new - n_avg) < tolerance:
            n_avg = n_avg_new
            break
        n_avg = n_avg_new

    # --- Output Results ---
    print("Mean-Field Self-Consistent Equation:")
    print("-------------------------------------")
    print("<n> = 1 / (1 + exp(-beta*(mu - z*epsilon*<n>)))")
    print("\nSubstituting the dimensionless parameters:")
    print(f"beta*mu = {beta_mu}")
    print(f"z = {z_horizontal} + {z_vertical} = {z}")
    print(f"-beta*epsilon*z = {beta_epsilon_z_term:.4f}")
    
    print("\nThe final numerical equation to be solved is:")
    print(f"<n> = 1 / (1 + exp(-({beta_mu} + {beta_epsilon_z_term:.4f} * <n>)))")

    print("\n--- Solution ---")
    print(f"The calculated average occupancy per site <n> is: {n_avg:.3f}")

solve_average_occupancy()
<<<0.848>>>