import math

def solve_average_occupancy():
    """
    Calculates the average occupancy per site using the mean-field approximation
    for a lattice gas model.
    """
    # --- Given Parameters ---
    z_horizontial = 4
    z_vertical = 8
    # mu = 0.1 * kB * T
    mu_over_kBT = 0.1
    # epsilon = -kB * T / (2 * pi)
    epsilon_over_kBT = -1 / (2 * math.pi)

    # --- Calculation ---
    # 1. Calculate the total coordination number
    z = z_horizontial + z_vertical

    # 2. Print the self-consistency equation with numerical values
    print("Under the mean-field approximation, the average occupancy <n> is found by solving the self-consistency equation:")
    print("<n> = 1 / (exp[(z * (epsilon/k_B*T) * <n> - (mu/k_B*T))] + 1)\n")
    print("Substituting the given parameter values:")
    print(f"Total coordination number z = {z_horizontial} + {z_vertical} = {z}")
    print(f"Interaction energy ratio epsilon/k_B*T = -1/(2*pi) ≈ {epsilon_over_kBT:.4f}")
    print(f"Chemical potential ratio mu/k_B*T = {mu_over_kBT}")
    
    # Using z * epsilon_over_kBT to show the full term in the equation
    z_eps_term = z * epsilon_over_kBT
    
    print(f"\nThe final equation to solve for <n> is:")
    print(f"<n> = 1 / (exp[{z_eps_term:.4f} * <n> - {mu_over_kBT}] + 1)\n")

    # 3. Solve the equation numerically using fixed-point iteration
    n_avg = 0.5  # Initial guess for <n>
    tolerance = 1e-9
    max_iterations = 100

    for i in range(max_iterations):
        n_avg_old = n_avg
        exponent = z_eps_term * n_avg - mu_over_kBT
        n_avg_new = 1.0 / (math.exp(exponent) + 1)
        
        # Check for convergence
        if abs(n_avg_new - n_avg_old) < tolerance:
            break
        n_avg = n_avg_new

    # 4. Print the final result
    print(f"After numerical iteration, the average occupancy <n> is: {n_avg:.3f}")
    
    return n_avg

# Execute the function and store the result
final_occupancy = solve_average_occupancy()