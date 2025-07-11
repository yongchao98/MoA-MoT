import numpy as np

def solve_occupancy():
    """
    Calculates the average occupancy per site <n> for a lattice gas system
    using the mean-field approximation and a numerical iterative solver.
    """
    # Define the parameters from the problem statement
    beta_mu = 0.1
    # Note: The expression for epsilon is ε = -k_B*T / (2*pi)
    # So, βε = (1/(k_B*T)) * (-k_B*T / (2*pi)) = -1/(2*pi)
    beta_epsilon = -1 / (2 * np.pi)
    z_h = 4
    z_v = 8
    z = z_h + z_v

    # The self-consistency equation is: <n> = 1 / (exp(z * βε * <n> - βμ) + 1)
    # We will print the equation with the calculated coefficients.
    print("Solving the self-consistency equation for the average occupancy <n>:")
    print("<n> = 1 / (exp(z * βε * <n> - βμ) + 1)\n")
    print("With the given parameters, the equation becomes:")
    
    coeff_n = z * beta_epsilon
    const_term = -beta_mu
    
    print(f"  z = {z}")
    print(f"  βε = {beta_epsilon:.4f}")
    print(f"  βμ = {beta_mu:.4f}\n")
    print("The final equation to solve numerically is:")
    print(f"<n> = 1 / (exp({coeff_n:.4f} * <n> {const_term:.4f}) + 1)\n")

    # Numerical solution using fixed-point iteration
    n = 0.5  # Initial guess for <n>
    tolerance = 1e-6
    max_iterations = 100
    
    for i in range(max_iterations):
        n_old = n
        exponent = z * beta_epsilon * n_old - beta_mu
        n = 1 / (np.exp(exponent) + 1)
        
        if np.abs(n - n_old) < tolerance:
            break

    # Output the final result
    print(f"The calculated average occupancy per site <n> is: {n:.3f}")
    return n

# Run the calculation and store the result
average_occupancy = solve_occupancy()
# The final answer must be returned in the specified format
final_answer = f"<<<{average_occupancy:.3f}>>>"
# print(final_answer) # This print is for the final answer block and should be commented out here
# but the format will be used at the very end of the response.