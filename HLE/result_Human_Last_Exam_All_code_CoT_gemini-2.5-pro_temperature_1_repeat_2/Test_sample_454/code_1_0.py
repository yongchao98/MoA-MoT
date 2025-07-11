import numpy as np
from scipy.optimize import fsolve

def solve_lattice_gas_occupancy():
    """
    Calculates the average occupancy per site for a lattice gas model
    using the mean-field approximation.
    """
    # 1. Define the parameters based on the problem statement.
    # The term k_B * T cancels out, so we can use dimensionless quantities.
    # beta_epsilon corresponds to βε
    beta_epsilon = -1 / (2 * np.pi)
    
    # beta_mu corresponds to βμ
    beta_mu = 0.1
    
    # The coordination number z is the sum of horizontal and vertical neighbors
    # for a particle in the bulk of the multilayer structure.
    z_horizontal = 4
    z_vertical = 8
    z = z_horizontal + z_vertical
    
    # 2. Define the self-consistency equation to be solved.
    # The equation is <n> = 1 / (exp(β(ε*z*<n> - μ)) + 1).
    # We rearrange it to a function f(<n>) = 0 for the solver.
    # f(<n>) = <n> - 1 / (exp(βε*z*<n> - βμ) + 1)
    def equation_to_solve(n):
        exponent = beta_epsilon * z * n - beta_mu
        return n - 1.0 / (np.exp(exponent) + 1.0)
        
    # 3. Solve the equation numerically.
    # We provide an initial guess for the solver. <n> is a probability, so 0.5 is a reasonable guess.
    initial_guess = 0.5
    # fsolve returns an array, so we take the first element.
    average_occupancy = fsolve(equation_to_solve, initial_guess)[0]
    
    # 4. Print the final results as requested.
    # The problem asks to output each number in the final equation.
    print("The self-consistency equation for the average occupancy <n> is:")
    print("<n> = 1 / (exp(β(ε*z*<n> - μ)) + 1)\n")
    print("Substituting the dimensionless parameters where β = 1/(k_B*T):")
    print("βε = -1 / (2*π) ≈ {:.4f}".format(beta_epsilon))
    print("βμ = 0.1")
    print("z = {} + {} = {}".format(z_horizontal, z_vertical, z))
    
    print("\nThe equation with numerical values is:")
    # Using format to have more control over the output string
    equation_str = "<n> = 1 / (exp({:.4f} * {} * <n> - {:.1f}) + 1)".format(beta_epsilon, z, beta_mu)
    print(equation_str)
    
    final_equation_str = "<n> = 1 / (exp({:.4f}*<n> - {:.1f}) + 1)".format(beta_epsilon * z, beta_mu)
    print(final_equation_str)
    
    print("\nSolving this transcendental equation numerically yields:")
    print("The average occupancy per site <n> is: {:.3f}".format(average_occupancy))

# Execute the function to get the result.
solve_lattice_gas_occupancy()
<<<0.848>>>