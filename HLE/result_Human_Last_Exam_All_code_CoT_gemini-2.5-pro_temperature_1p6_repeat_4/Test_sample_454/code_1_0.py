import math

def solve_mean_field_occupancy():
    """
    Calculates the average occupancy per site for a lattice gas model
    using the mean-field approximation.
    """
    # Given parameters
    mu_over_kBT = 0.1
    # epsilon_over_kBT is not explicitly a variable, but is part of the coefficient
    z_horizontal = 4
    z_vertical = 8
    
    # Total coordination number for a site in a non-monolayer
    z = z_horizontal + z_vertical
    
    # The self-consistency equation is of the form:
    # <n> = 1 / (1 + exp(z * beta*epsilon * <n> - beta*mu))
    # where beta*epsilon = -1 / (2*pi) and beta*mu = 0.1
    # <n> = 1 / (1 + exp(-z/(2*pi) * <n> - 0.1))
    
    print("Solving the self-consistency equation for the average occupancy <n>:")
    print(f"<n> = 1 / (1 + exp(-({z} / (2 * pi)) * <n> - {mu_over_kBT}))")
    print("-" * 20)

    # Initial guess for the average occupancy <n>
    n_avg = 0.5 
    
    # Coefficient for <n> in the exponent
    coeff = -z / (2 * math.pi)
    
    # Iterate to find the solution for <n>
    # 50 iterations are more than sufficient for convergence
    for _ in range(50):
        exponent = coeff * n_avg - mu_over_kBT
        n_avg = 1 / (1 + math.exp(exponent))
        
    print(f"The calculated average occupancy per site <n> is: {n_avg:.3f}")

solve_mean_field_occupancy()
<<<0.848>>>