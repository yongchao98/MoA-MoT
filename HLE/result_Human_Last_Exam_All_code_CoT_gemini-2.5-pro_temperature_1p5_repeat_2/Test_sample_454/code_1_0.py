import math

def solve_lattice_gas():
    """
    Calculates the average occupancy per site for a lattice gas model
    using the mean-field approximation and fixed-point iteration.
    """
    # Define parameters based on the problem statement
    # Chemical potential over kT
    mu_over_kT = 0.1
    # Interaction energy over kT
    eps_over_kT = -1.0 / (2.0 * math.pi)
    # Coordination numbers
    z_horizontal = 4
    z_vertical = 8
    # Total coordination number for the mean-field approximation
    z = z_horizontal + z_vertical

    # The self-consistency equation to solve is:
    # <n> = 1 / (1 + exp((z * eps/kT * <n> - mu/kT)))
    # We use a fixed-point iteration method.

    # Initial guess for <n> (must be between 0 and 1)
    n_avg = 0.5

    # Iterate to find the solution (100 iterations is sufficient for convergence)
    for _ in range(100):
        exponent = z * eps_over_kT * n_avg - mu_over_kT
        n_avg = 1.0 / (1.0 + math.exp(exponent))

    # The final converged value of <n>
    final_n_avg = n_avg

    # Round the values for the final output equation
    final_n_avg_rounded = round(final_n_avg, 3)
    eps_over_kT_rounded = round(eps_over_kT, 3)
    mu_over_kT_rounded = round(mu_over_kT, 3)
    z_val = int(z)

    # Print the explanation and the final solved equation with numerical values
    print("The self-consistency equation for the average occupancy <n> is solved:")
    print("<n> = 1 / (1 + exp((z * (epsilon/kT) * <n> - (mu/kT))))\n")
    print("Using the given parameters:")
    print(f"z = {z_val}")
    print(f"epsilon/kT = {eps_over_kT_rounded}")
    print(f"mu/kT = {mu_over_kT_rounded}\n")
    print(f"The numerically calculated average occupancy <n> is: {final_n_avg_rounded}\n")
    print("The final equation with the solution <n> plugged in is:")
    # The f-string formats the final equation with all numbers included
    print(f"{final_n_avg_rounded} = 1 / (1 + exp(({z_val} * {eps_over_kT_rounded} * {final_n_avg_rounded}) - {mu_over_kT_rounded}))")

if __name__ == '__main__':
    solve_lattice_gas()