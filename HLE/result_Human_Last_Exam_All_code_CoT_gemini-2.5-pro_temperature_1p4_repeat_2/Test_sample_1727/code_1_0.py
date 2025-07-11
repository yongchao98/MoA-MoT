import numpy as np
from scipy.optimize import fsolve

def calculate_chi_at_max_response():
    """
    Calculates the stable amplitude factor chi by finding the maximum
    response of the system.
    
    The steps are:
    1. Define the equation whose root gives the maximum chi. The maximum
       chi occurs when H satisfies tanh(H) = 1/H.
    2. Use a numerical solver to find the root H_max.
    3. Substitute H_max into the simplified expression for chi at the maximum.
    4. Print the final result.
    """
    
    # 1. Define the function for the root-finding problem: f(H) = tanh(H) - 1/H
    # We want to find H such that f(H) = 0.
    equation_to_solve = lambda H: np.tanh(H) - 1/H
    
    # 2. Use a numerical solver (fsolve) to find the root.
    # An initial guess of 1.0 is reasonable, as tanh(1) ~ 0.76 and 1/1=1.
    initial_guess = 1.0
    H_max_solution = fsolve(equation_to_solve, initial_guess)
    H_max = H_max_solution[0]

    # 3. Calculate chi using the simplified formula for the maximum:
    # chi_max = sqrt(H_max^2 - 1) / H_max^2
    chi = np.sqrt(H_max**2 - 1) / H_max**2
    
    # The problem defines the relationship e_s = chi * S0 / (delta_rho * g)
    # We are providing the calculated value for the factor chi.
    print(f"The dimensionless parameter H=kh that maximizes the response is: {H_max:.4f}")
    print(f"The stable amplitude factor chi is: {chi:.4f}")

    # As requested, outputting each number in the final equation.
    # The final equation describes the surface topography e_s based on our calculated chi.
    print("\nThe equation for the stable vertical displacement (e_s) is:")
    e_s_val = chi
    S0_val = 1.0  # representing S0 as a variable
    delta_rho_g_val = 1.0 # representing the product of delta_rho and g
    print(f"e_s = {e_s_val:.4f} * (S0 / (delta_rho * g))")

calculate_chi_at_max_response()