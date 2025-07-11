import numpy as np
from scipy.optimize import fsolve

def solve_adsorption():
    """
    Solves the self-consistent mean-field equation for monolayer adsorption.
    """
    # Parameters given in the problem (in units of k_B*T)
    epsilon_1 = 0.1
    mu = 0.15
    z_l = 4
    
    # We interpret the problem for k=1 (monolayer)
    k = 1
    epsilon_l = (0.02)**k
    
    # The self-consistent equation for theta is:
    # ln(theta / (1 - theta)) = mu + epsilon_1 + z_l * epsilon_l * theta
    # f(theta) = ln(theta / (1 - theta)) - (z_l * epsilon_l * theta) - (mu + epsilon_1) = 0
    
    # Define the function to find the root of
    def equation(theta):
        if theta <= 0 or theta >= 1:
            return 1e9 # Return a large number outside the valid domain
        return np.log(theta / (1 - theta)) - z_l * epsilon_l * theta - (mu + epsilon_1)

    # Initial guess for theta
    initial_guess = 0.5
    
    # Solve the equation numerically
    theta_solution = fsolve(equation, initial_guess)[0]
    
    # The average number of adsorbed layers per site is theta for a monolayer
    avg_layers = theta_solution
    
    print("Based on the interpretation that the problem parameters imply a monolayer system (k=1):")
    print("The self-consistent equation to solve for the coverage theta is:")
    print(f"ln(theta / (1 - theta)) = mu + epsilon_1 + z_l * epsilon_l * theta")
    print("Plugging in the numeric values (in units of k_B*T):")
    print(f"ln(theta / (1 - theta)) = {mu} + {epsilon_1} + {z_l} * {epsilon_l} * theta")
    final_eq = f"ln(theta / (1 - theta)) - {z_l * epsilon_l:.2f}*theta - {mu + epsilon_1:.2f} = 0"
    print(f"Final equation to solve: {final_eq}")
    print("\nSolving this equation numerically gives the average number of adsorbed layers per site:")
    print(f"{avg_layers}")
    
solve_adsorption()