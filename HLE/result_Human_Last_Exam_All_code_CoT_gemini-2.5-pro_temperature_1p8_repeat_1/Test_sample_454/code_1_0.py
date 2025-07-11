import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site using the mean-field approximation.
    """
    # Parameters given in the problem
    # epsilon / (k_B * T)
    beta_epsilon = -1.0 / (2.0 * math.pi)
    # mu / (k_B * T)
    beta_mu = 0.1
    z_horizontal = 4
    z_vertical = 8
    
    # Total coordination number for a particle in the bulk
    z = z_horizontal + z_vertical
    
    # Self-consistent equation for average occupancy <n>:
    # <n> = 1 / (exp(z * beta_epsilon * <n> - beta_mu) + 1)
    
    # We solve this using fixed-point iteration: n_k+1 = f(n_k)
    # Initial guess for <n>
    n = 0.5
    
    # Iteration settings
    max_iterations = 100
    tolerance = 1e-7
    
    for i in range(max_iterations):
        exponent = z * beta_epsilon * n - beta_mu
        n_new = 1.0 / (math.exp(exponent) + 1.0)
        
        if abs(n_new - n) < tolerance:
            break
        n = n_new

    # Output the components of the equation
    print("This code solves the self-consistent mean-field equation for the average occupancy <n>:")
    print("<n> = 1 / (exp(z * (epsilon/kT) * <n> - (mu/kT)) + 1)\n")
    print("Using the following parameters:")
    print(f"z = z_horizontal + z_vertical = {z_horizontal} + {z_vertical} = {z}")
    print(f"epsilon/kT = -1/(2*pi) approx {beta_epsilon:.5f}")
    print(f"mu/kT = {beta_mu:.5f}\n")
    
    # Output the final equation with numerical values
    final_equation = f"<n> = 1 / (exp({z} * ({beta_epsilon:.5f}) * <n> - {beta_mu:.5f}) + 1)"
    print("The final equation being solved is:")
    print(final_equation)
    
    # Output the result
    print(f"\nThe calculated average occupancy per site <n> is: {n:.3f}")

solve_average_occupancy()