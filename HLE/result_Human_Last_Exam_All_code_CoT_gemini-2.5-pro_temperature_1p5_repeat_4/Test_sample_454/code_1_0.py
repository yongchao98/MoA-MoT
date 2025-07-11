import math

def calculate_average_occupancy():
    """
    Calculates the average site occupancy <n> for a lattice gas model
    using the mean-field approximation and solves the resulting
    self-consistency equation numerically.
    """
    # Parameters from the problem description
    z_horizontial = 4
    z_vertical = 8
    
    # Dimensionless parameters (energies in units of k_B*T)
    # mu' = mu / (k_B * T)
    beta_mu = 0.1
    # epsilon' = epsilon / (k_B * T)
    beta_eps = -1.0 / (2.0 * math.pi)

    # The total effective coordination number
    z = z_horizontial + z_vertical

    # --- Numerical Solution using Fixed-Point Iteration ---
    # The self-consistency equation is: n_avg = f(n_avg)
    # where f(x) = 1 / (exp(z * beta_eps * x - beta_mu) + 1)
    
    n_avg = 0.5  # Initial guess for the average occupancy
    tolerance = 1e-9 # Convergence criterion

    for _ in range(100):  # Iterate up to 100 times to ensure convergence
        exponent = z * beta_eps * n_avg - beta_mu
        n_avg_new = 1.0 / (math.exp(exponent) + 1.0)
        
        # Check if the value has converged
        if abs(n_avg_new - n_avg) < tolerance:
            n_avg = n_avg_new
            break
        n_avg = n_avg_new

    # --- Output Results as specified ---
    print("This script solves for the average site occupancy <n> in a mean-field lattice gas model.")
    print("The self-consistency equation is: <n> = 1 / (exp(z * (eps/kBT) * <n> - mu/kBT) + 1)\n")
    
    print("The numerical values for each parameter in the equation are:")
    print(f"Total coordination number, z = {z}")
    print(f"Dimensionless interaction energy, eps/kBT = {beta_eps:.4f}")
    print(f"Dimensionless chemical potential, mu/kBT = {beta_mu:.1f}\n")

    # Output the final equation with the calculated value of <n> plugged in
    print("The final equation with the calculated numbers is:")
    print(f"{n_avg:.3f} ≈ 1 / (exp({z} * ({beta_eps:.4f}) * {n_avg:.3f} - {beta_mu:.1f}) + 1)")

    # Output the final numerical answer
    print(f"\nThe calculated average occupancy per site, <n>, is: {n_avg:.3f}")


if __name__ == '__main__':
    calculate_average_occupancy()