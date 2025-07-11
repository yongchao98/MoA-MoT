import math

def solve_lattice_gas():
    """
    Calculates the average occupancy per site for a lattice gas system
    using the mean-field approximation.
    """
    # 1. Define the given parameters.
    # The problem can be solved using dimensionless ratios, so we don't need k_B or T explicitly.
    z_horizontal = 4
    z_vertical = 8
    mu_over_kBT = 0.1
    epsilon_over_kBT = -1 / (2 * math.pi)

    # The total coordination number z for a site in the bulk of a multilayer system
    # is the sum of its in-plane (horizontal) and out-of-plane (vertical) neighbors.
    z = z_horizontal + z_vertical

    # 2. Solve the self-consistency equation numerically.
    # The equation is: <n> = 1 / (exp(beta * (z * epsilon * <n> - mu)) + 1)
    # Let x = <n>, the average occupancy. The equation becomes:
    # x = 1 / (exp(z * (epsilon/kBT) * x - (mu/kBT)) + 1)

    x = 0.5  # Initial guess for the average occupancy <n>

    # We iterate to find the fixed point, which is the solution for x.
    # 100 iterations are more than enough for convergence here.
    for _ in range(100):
        exponent = z * epsilon_over_kBT * x - mu_over_kBT
        x = 1 / (math.exp(exponent) + 1)

    # 3. Round the final result to three decimal places.
    average_occupancy = round(x, 3)

    # 4. Print the explanation, the final equation with numbers, and the result.
    print("The self-consistency equation for the average occupancy <n> is:")
    print("<n> = 1 / (exp(z * (epsilon / (k_B*T)) * <n> - mu / (k_B*T)) + 1)\n")
    
    print("Using the provided values, the numbers in the final equation are:")
    print(f"Total coordination number, z = {z}")
    print(f"Dimensionless interaction energy, epsilon / (k_B*T) = {epsilon_over_kBT:.3f}")
    print(f"Dimensionless chemical potential, mu / (k_B*T) = {mu_over_kBT:.3f}\n")

    print("The final equation being solved numerically is:")
    z_times_eps_eff = round(z * epsilon_over_kBT, 3)
    print(f"<n> = 1 / (exp({z_times_eps_eff} * <n> - {mu_over_kBT}) + 1)\n")
    
    print(f"The calculated average occupancy per site is:")
    print(f"<n> = {average_occupancy}")

    print(f"\n<<<{average_occupancy}>>>")

# Execute the function
solve_lattice_gas()