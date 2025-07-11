import numpy as np

def calculate_force_expression():
    """
    This function prints the components and the final expression for the force
    per unit area based on the selected answer choice.
    """
    # Define the symbolic components of the equation
    direction = "i_x"
    factor_half = "1/2"
    mu_0 = "mu_0"
    K_0_sq = "K_0^2"
    time_term = "cos^2(omega * t)"
    cosh_term_denom = "cosh^2((omega_p * d) / c)"
    exp_term = "exp(-(omega * d) / c)"

    # Construct the final expression string
    final_expression = (
        f"f_vector = {direction} * {factor_half} * "
        f"({mu_0} * {K_0_sq} * {time_term} / {cosh_term_denom}) * {exp_term}"
    )

    print("The final expression for the force per unit area is constructed from the following components:")
    print(f"Direction vector: {direction}")
    print(f"Numerical factor: {factor_half}")
    print(f"Permeability of free space: {mu_0}")
    print(f"Square of surface current amplitude: {K_0_sq}")
    print(f"Time-dependent term: {time_term}")
    print(f"Denominator term from field decay: {cosh_term_denom}")
    print(f"Additional exponential factor: {exp_term}\n")
    
    print("Final equation:")
    # We print each part of the equation as requested
    print(f"   f = (i_x) * (1/2) * mu_0 * K_0^2 * (cos^2(omega*t) / cosh^2((omega_p*d)/c)) * exp(-(omega*d)/c)")


# Execute the function to display the result
calculate_force_expression()
