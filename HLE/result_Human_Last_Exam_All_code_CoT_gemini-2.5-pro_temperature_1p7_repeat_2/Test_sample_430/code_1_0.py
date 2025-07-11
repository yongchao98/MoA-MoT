import numpy as np
from scipy.optimize import brentq

def solve_adsorption_model():
    """
    Solves the lattice adsorption problem under the assumption of monolayer adsorption (k=1).
    """
    # Assumption: The problem describes monolayer adsorption (k=1).
    # This assumption is made because crucial parameters for multilayer adsorption 
    # (e.g., epsilon_2, epsilon_inter) are not provided, making the multilayer 
    # case unsolvable as stated.
    k = 1

    # Parameters from the problem statement (in dimensionless units of beta = 1/(kB*T))
    beta_epsilon_1 = 0.1
    beta_mu = 0.15
    z_l = 4
    # Calculate beta*epsilon_l based on the assumption k=1
    beta_epsilon_l = (0.02)**k
    
    # Define the function f(theta) = 0 that we need to solve.
    # This is the rearranged Fowler-Guggenheim self-consistency equation:
    # ln(theta / (1 - theta)) - (beta*mu + beta*epsilon_1 + z_l * beta*epsilon_l * theta) = 0
    def equation_to_solve(theta, mu, eps1, zl, epsl):
        # We need to solve for theta in [0, 1].
        # The log term is handled by the solver's range.
        log_term = np.log(theta / (1.0 - theta))
        energy_term = mu + eps1
        interaction_term = zl * epsl * theta
        return log_term - energy_term - interaction_term

    # Print the explanation and the specific equation being solved
    print("Based on the assumption of a monolayer model (k=1), the system is described by the Fowler-Guggenheim isotherm.")
    print("The self-consistency equation for the fractional coverage theta is:")
    print("ln(theta / (1 - theta)) = beta*mu + beta*epsilon_1 + z_l * beta*epsilon_l * theta\n")

    print("Plugging in the numerical values provided:")
    print(f"ln(theta / (1 - theta)) = {beta_mu} + {beta_epsilon_1} + {z_l} * {beta_epsilon_l} * theta")
    
    # Consolidate the constant terms for a cleaner final equation
    constant_term = beta_mu + beta_epsilon_1
    theta_coeff = z_l * beta_epsilon_l
    print(f"ln(theta / (1 - theta)) - {theta_coeff:.2f}*theta - {constant_term:.2f} = 0\n")

    # Numerically solve for theta.
    # The coverage theta must be between 0 and 1. We use a small tolerance 
    # to keep the argument of the logarithm positive.
    try:
        solution_theta = brentq(
            equation_to_solve, 
            a=1e-9, 
            b=1-1e-9, 
            args=(beta_mu, beta_epsilon_1, z_l, beta_epsilon_l)
        )
        print("Solving this equation for theta gives the average number of adsorbed layers per site.")
        print(f"The result is: {solution_theta}")
        # Final answer format
        print(f"<<< {solution_theta} >>>")
    except ValueError:
        print("A solution for theta could not be found in the interval (0, 1).")

# Execute the solution
solve_adsorption_model()