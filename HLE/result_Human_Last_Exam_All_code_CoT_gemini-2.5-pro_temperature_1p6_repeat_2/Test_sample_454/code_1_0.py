import math

def solve_mean_field_occupancy():
    """
    Solves the self-consistency equation for average lattice site occupancy
    using the mean-field approximation.
    """
    # 1. Define the parameters based on the problem statement
    # We work with dimensionless quantities epsilon' = epsilon / (k_B*T) and mu' = mu / (k_B*T)
    epsilon_prime = -1.0 / (2.0 * math.pi)
    mu_prime = 0.1
    z_horizontial = 4
    z_vertical = 8
    z_total = z_horizontial + z_vertical

    # 2. Define constants for the self-consistency equation: <n> = 1 / (exp(A*<n> + B) + 1)
    # A = z * epsilon / (k_B*T)
    # B = -mu / (k_B*T)
    A = z_total * epsilon_prime
    B = -mu_prime
    
    # 3. Print the derived equation with numerical coefficients
    print("Based on the mean-field approximation, we solve the self-consistency equation for the average occupancy <n>:")
    print("<n> = 1 / (exp(A * <n> + B) + 1)\n")
    print(f"Using the given parameters, the coefficients are:")
    print(f"Total coordination number z = {z_total}")
    # The final equation is <n> = 1/(exp((z*epsilon/kBT)*<n> - mu/kBT)+1)
    print(f"z * epsilon / (k_B * T) = {z_total} * ({epsilon_prime:.4f}) = {A:.4f}")
    print(f"-mu / (k_B * T) = {-mu_prime:.4f}")
    
    # 4. Solve the equation numerically using an iterative method
    n = 0.5  # Initial guess for <n>
    tolerance = 1e-7
    for _ in range(100):  # Iterate to find a stable solution
        n_new = 1.0 / (math.exp(A * n + B) + 1.0)
        if abs(n_new - n) < tolerance:
            break
        n = n_new

    # 5. Print the final result
    print(f"\nSolving the equation {n:.3f} = 1 / (exp({A:.3f} * {n:.3f} + {B:.3f}) + 1) numerically gives:")
    print(f"The average occupancy per site <n> is: {n:.3f}")

# Execute the function to find the solution
solve_mean_field_occupancy()
<<<0.848>>>