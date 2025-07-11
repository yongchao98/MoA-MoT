import numpy as np
from scipy.optimize import fsolve

def solve_occupancy():
    """
    Solves for the average occupancy <n> based on the mean-field model
    for multilayer adsorption, assuming a repulsive interaction energy to ensure
    a physical solution.
    """

    # The self-consistent equation to solve is of the form:
    # mu_eff - z * epsilon_eff * y = ln(y / (1 + y))
    # where y = <n>, mu_eff = mu / (k_B*T) = 0.1, and z*epsilon_eff = 12 * (1/(2*pi)) = 6/pi.
    
    # We define a function whose root we want to find, f(y) = 0
    # f(y) = mu_eff - z*epsilon_eff*y - ln(y/(1+y))
    def equation(y, mu_eff, z_epsilon_eff):
        if y <= 0:
            # The function is undefined for y<=0, return a large value
            # to guide the solver away from this region.
            return 1e9
        return mu_eff - z_epsilon_eff * y - np.log(y / (1 + y))

    # Dimensionless parameters from the problem description
    # We assume epsilon is positive, so the interaction is repulsive.
    mu_eff = 0.1
    z_eff = 12.0
    epsilon_eff = 1.0 / (2.0 * np.pi)
    z_epsilon_eff = z_eff * epsilon_eff
    
    # Initial guess for the root <n>
    initial_guess = 0.5
    
    # Solve the equation numerically
    n_average = fsolve(equation, initial_guess, args=(mu_eff, z_epsilon_eff))[0]

    # Print the explanation and the final equation being solved
    print("The derived self-consistent equation for the average occupancy <n> is:")
    print("β(μ - zε<n>) = ln(<n> / (1 + <n>))")
    print("\nUsing the given parameters results in a contradiction (positive = negative).")
    print("Assuming a typo in the sign of the interaction energy ε to enable a solution.")
    print("We proceed with ε = +k_B*T/(2π) instead of -k_B*T/(2π).")
    print("\nThe equation in dimensionless form (dividing by k_B*T) is:")
    print(f"{mu_eff:.3f} - ({z_eff:.3f} * {epsilon_eff:.3f}) * <n> = ln(<n> / (1 + <n>))")
    print(f"{mu_eff:.3f} - {z_epsilon_eff:.3f} * <n> = ln(<n> / (1 + <n>))")
    
    # Print the final result
    print(f"\nThe numerically solved average occupancy per site <n> is: {n_average:.3f}")
    
solve_occupancy()
<<<0.578>>>