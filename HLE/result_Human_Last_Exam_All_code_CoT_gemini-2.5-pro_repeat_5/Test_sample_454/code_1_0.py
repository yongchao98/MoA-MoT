import numpy as np
from scipy.optimize import fsolve

def solve_lattice_gas():
    """
    Calculates the average occupancy per site for a lattice gas
    using the mean-field approximation.
    """
    # Define dimensionless parameters based on the problem statement
    # mu = 0.1 * k_B * T  => beta * mu = 0.1
    beta_mu = 0.1

    # epsilon = -k_B * T / (2 * pi) => beta * epsilon = -1 / (2 * pi)
    beta_epsilon = -1 / (2 * np.pi)

    # Effective coordination number z = z_horizontal + z_vertical
    z_h = 4
    z_v = 8
    z = z_h + z_v

    # The self-consistency equation is of the form:
    # n = 1 / (exp(-(beta*mu - z*beta*epsilon*n)) + 1)
    # We define a function f(n) = 0 to solve it numerically.
    def self_consistency_function(n):
        # n represents the average occupancy <n>
        exponent = -(beta_mu - z * beta_epsilon * n)
        return n - 1 / (np.exp(exponent) + 1)

    # Provide an initial guess for the solver (occupancy is between 0 and 1)
    initial_guess = 0.5

    # Use fsolve to find the root of the function
    solution = fsolve(self_consistency_function, initial_guess)
    avg_occupancy = solution[0]

    # Print the final equation with the numerical values plugged in
    print("The self-consistency equation for the average occupancy <n> is:")
    print(f"<n> = 1 / (exp(-({beta_mu} - {z} * ({beta_epsilon:.4f}) * <n>)) + 1)")
    print(f"<n> = 1 / (exp(-({beta_mu} + {abs(z * beta_epsilon):.4f} * <n>)) + 1)")
    
    # Print the final result
    print("\nSolving this equation numerically gives the result:")
    print(f"The average occupancy per site <n> is: {avg_occupancy:.3f}")

if __name__ == "__main__":
    solve_lattice_gas()