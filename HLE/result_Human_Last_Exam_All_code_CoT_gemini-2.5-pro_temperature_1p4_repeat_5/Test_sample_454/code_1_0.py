import math

def solve_occupancy():
    """
    Calculates the average occupancy per site using the mean-field approximation.
    """
    # Parameters from the problem description
    mu_over_kT = 0.1
    epsilon_over_kT = -1 / (2 * math.pi)
    z_horizontial = 4
    z_vertical = 8
    
    # In mean-field theory, we use the total coordination number
    z = z_horizontial + z_vertical
    
    # The term beta * z * epsilon
    beta_z_epsilon = z * epsilon_over_kT
    
    # The term beta * mu
    beta_mu = mu_over_kT

    # Print the explanation and the equation with numerical values
    print("The self-consistent mean-field equation for the average occupancy ⟨n⟩ is:")
    print("⟨n⟩ = 1 / (exp(β(εz⟨n⟩ - μ)) + 1)")
    print("\nWith the given parameters, the terms in the exponent are:")
    print(f"βμ = {beta_mu:.3f}")
    print(f"βzε = {beta_z_epsilon:.3f}")
    print("\nThus, the equation to solve for ⟨n⟩ is:")
    print(f"⟨n⟩ = 1 / (exp({beta_z_epsilon:.3f} * ⟨n⟩ - {beta_mu:.3f}) + 1)")

    # Numerically solve the transcendental equation using fixed-point iteration
    # Let x = ⟨n⟩, the equation is x = f(x)
    # Start with an initial guess for ⟨n⟩
    n_avg = 0.5 
    
    # Iterate to find the solution
    for _ in range(100): # 100 iterations are more than enough for convergence
        exponent = beta_z_epsilon * n_avg - beta_mu
        n_avg_new = 1 / (math.exp(exponent) + 1)
        # Check for convergence
        if abs(n_avg_new - n_avg) < 1e-6:
            break
        n_avg = n_avg_new
        
    # Print the final result
    print("\nSolving this equation numerically yields the average occupancy per site:")
    print(f"⟨n⟩ = {n_avg:.3f}")

if __name__ == "__main__":
    solve_occupancy()