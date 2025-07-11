import math

def solve_mean_field_occupancy():
    """
    Calculates the average occupancy per site <n> for a lattice gas
    using the mean-field approximation.
    """
    # Parameters given in the problem
    # mu = 0.1 * k_B * T  => beta * mu = 0.1
    # epsilon = -k_B * T / (2 * pi) => beta * epsilon = -1 / (2 * pi)
    # z_horizontial = 4 (coordination number for the 2D layer)
    
    beta_mu = 0.1
    z = 4.0
    # Note: pi is imported from the math module
    beta_epsilon = -1.0 / (2.0 * math.pi)

    # The self-consistency equation is:
    # <n> = 1 / (1 + exp[-(beta*mu - z * beta*epsilon * <n>)])
    # Let x = <n>
    # x = 1 / (1 + exp[-(0.1 - 4 * (-1/(2*pi)) * x)])
    # x = 1 / (1 + exp[-(0.1 + (2/pi) * x)])
    
    print("Solving the mean-field self-consistency equation for the average occupancy <n>:")
    print(f"<n> = 1 / (1 + exp[ -({beta_mu} - ({z}) * ({beta_epsilon:.4f}) * <n>) ])")
    print("which simplifies to:")
    print(f"<n> = 1 / (1 + exp[ -({beta_mu} + {(-z*beta_epsilon):.4f} * <n>) ])\n")


    # Iterative solver for the transcendental equation x = f(x)
    
    # Initial guess for <n>
    n = 0.5 
    
    # Tolerance for convergence
    tolerance = 1e-7
    
    # Maximum number of iterations to prevent infinite loops
    max_iterations = 100
    
    for i in range(max_iterations):
        # The right-hand side of the equation
        n_next = 1.0 / (1.0 + math.exp(-(beta_mu - z * beta_epsilon * n)))
        
        # Check for convergence
        if abs(n_next - n) < tolerance:
            n = n_next
            break
        
        n = n_next
    else:
        # This part is executed if the loop finishes without break
        print("Warning: Solver did not converge within the maximum number of iterations.")

    # Print the final result rounded to three decimal places
    print(f"The calculated average occupancy per site <n> is: {n:.3f}")

    return n

if __name__ == '__main__':
    average_occupancy = solve_mean_field_occupancy()
    print(f"\n<<< {average_occupancy:.3f} >>>")
