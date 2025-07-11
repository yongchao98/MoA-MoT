import numpy as np
from scipy.optimize import fsolve

def solve_adsorption_occupancy():
    """
    Solves for the average occupancy in a lattice gas model using the
    mean-field approximation.
    """
    # Define dimensionless parameters based on the problem statement
    # mu / (k_B * T)
    mu_over_kT = 0.1
    # epsilon / (k_B * T)
    eps_over_kT = -1 / (2 * np.pi)
    
    # Coordination numbers
    z_h = 4
    z_v = 8
    
    # Based on the problem description, we interpret the system as a multilayer
    # adsorbate and calculate the occupancy for a "bulk" layer.
    # The total coordination number z is the sum of horizontal and vertical neighbors.
    z = z_h + z_v
    
    # Define the self-consistency equation to be solved.
    # We are looking for the root of the function f(x) - x = 0, where x is <n>.
    def self_consistency_equation(x, mu_kt, z_val, eps_kt):
        """
        Represents the equation: RHS - x = 0, where RHS is the mean-field expression for <n>.
        """
        exponent = -(mu_kt - z_val * eps_kt * x)
        rhs = 1.0 / (1.0 + np.exp(exponent))
        return rhs - x

    # We use a numerical solver to find the root of the equation.
    # An initial guess for the occupancy (between 0 and 1) is required.
    initial_guess = 0.5
    average_occupancy, = fsolve(self_consistency_equation, initial_guess, args=(mu_over_kT, z, eps_over_kT))

    # Output the steps and the result
    print("Solving the mean-field self-consistency equation for the average occupancy <n>:")
    print("<n> = 1 / (1 + exp( - (mu/(kT) - z * epsilon/(kT) * <n>) ))")
    print("\nGiven parameters lead to the following dimensionless values:")
    print(f"mu/(kT) = {mu_over_kT}")
    print(f"epsilon/(kT) = {eps_over_kT:.4f}")
    print(f"z_horizontial = {z_h}")
    print(f"z_vertical = {z_v}")
    print(f"Total coordination number used, z = {z_h} + {z_v} = {z}")
    
    print("\nThe specific equation to solve is:")
    # Note: z * eps_over_kT is negative, so the sign inside the bracket flips
    print(f"<n> = 1 / (1 + exp( -({mu_over_kT} + {abs(z * eps_over_kT):.4f} * <n>) ))")

    print(f"\nSolving for <n> numerically...")
    print(f"\nThe average occupancy per site <n> is: {average_occupancy:.3f}")

solve_adsorption_occupancy()
<<<0.848>>>