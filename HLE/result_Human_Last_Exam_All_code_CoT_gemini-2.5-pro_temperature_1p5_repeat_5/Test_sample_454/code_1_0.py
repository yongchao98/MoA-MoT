import math

def calculate_average_occupancy():
    """
    Calculates the average occupancy per site using the mean-field approximation.
    """
    # 1. Define the parameters based on the problem statement.
    z_horizontal = 4
    z_vertical = 8
    
    # In the mean-field model, a particle interacts with an average environment.
    # We assume the total coordination number 'z' is the sum of horizontal and
    # vertical neighbors for a representative particle in the adsorbed phase.
    z = z_horizontal + z_vertical

    # 2. Calculate the dimensionless energy and chemical potential terms.
    # The term beta = 1/(k_B*T) is implicit in these calculations.
    # beta_epsilon = beta * epsilon = (1/(k_B*T)) * (-k_B*T / (2*pi))
    beta_epsilon = -1.0 / (2.0 * math.pi)
    
    # beta_mu = beta * mu = (1/(k_B*T)) * (0.1 * k_B*T)
    beta_mu = 0.1

    # 3. Solve the self-consistency equation numerically.
    # The equation is: n = 1 / (1 + exp(z * beta_epsilon * n - beta_mu))
    # We use fixed-point iteration to find the value of n.
    n_avg = 0.5  # Initial guess for the average occupancy <n>
    
    # Iterate a sufficient number of times for the value to converge.
    for _ in range(100):
        exponent = z * beta_epsilon * n_avg - beta_mu
        n_avg = 1.0 / (1.0 + math.exp(exponent))

    # 4. Print the components of the final equation and the result.
    print("The final equation being solved is derived from the mean-field approximation:")
    print("<n> = 1 / (1 + exp(z * (β*ε) * <n> - (β*μ)))")
    
    print("\nSubstituting the given parameters:")
    print(f"Total coordination number z = {z_horizontal} + {z_vertical} = {z}")
    print(f"Dimensionless interaction term (β*ε) = -1 / (2*π) ≈ {beta_epsilon:.4f}")
    print(f"Dimensionless chemical potential term (β*μ) = {beta_mu}")
    
    print("\nFinal numerical equation solved:")
    print(f"<n> = 1 / (1 + exp({z} * {beta_epsilon:.4f} * <n> - {beta_mu}))")
    
    # Output the final calculated average occupancy rounded to three decimal places.
    print(f"\nThe average occupancy per site <n> is: {n_avg:.3f}")

calculate_average_occupancy()