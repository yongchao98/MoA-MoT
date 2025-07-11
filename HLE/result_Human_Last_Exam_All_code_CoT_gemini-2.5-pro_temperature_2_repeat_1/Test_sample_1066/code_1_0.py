import numpy as np
from scipy.special import digamma

def calculate_lower_bound_coefficients():
    """
    Calculates the coefficients for the lower bound on E[S].
    The bound has the form: E[S] >= n * alpha + n * C
    This function computes C.
    """
    # The Euler-Mascheroni constant gamma
    # Can be calculated as -digamma(1)
    gamma = -digamma(1)

    # The constant pi
    pi = np.pi

    # The constant from the Basel problem, zeta(2)
    pi_squared_over_6 = (pi**2) / 6

    # The coefficient for the n * alpha term is 1.0
    alpha_coeff = 1.0

    # The constant C in the n*C term of the lower bound
    constant_term_coeff = gamma - pi_squared_over_6

    print("The derived lower bound for E[S] is of the form:")
    print("E[S] >= A * n * alpha + B * n")
    print("\nCalculated coefficients:")
    print(f"A = {alpha_coeff}")
    print(f"B = gamma - pi^2/6 = {gamma} - {pi_squared_over_6} = {constant_term_coeff}")
    print("\nThus, the final inequality is:")
    print(f"E[S] >= {alpha_coeff} * n * alpha + {constant_term_coeff} * n")
    
    # Returning the final numerical value for the constant for the prompt format
    return constant_term_coeff

# Execute the function to find and print the result.
# We extract the final numerical coefficient B for the requested output format.
result_coefficient = calculate_lower_bound_coefficients()