import math

def calculate_optimal_steps_for_two_step_gd():
    """
    This function calculates and prints the optimal step sizes for a 
    two-step gradient descent method based on the Chebyshev polynomial approach.
    It demonstrates the calculation for a specific condition number M.
    """
    # Per the problem, we use M as the condition number and assume mu=1.
    # Let's use an example value for M.
    M = 10.0
    mu = 1.0

    print(f"Calculating the optimal step sizes for M = {M} and mu = {mu}:")
    print("-" * 50)

    # The formulas for the optimal step sizes are derived from the roots of a
    # shifted Chebyshev polynomial. The general formulas are:
    # gamma_1 = 2 / ( (M + mu) - (M - mu) / sqrt(2) )
    # gamma_2 = 2 / ( (M + mu) + (M - mu) / sqrt(2) )
    # The order of gamma_1 and gamma_2 is interchangeable.

    sqrt2_val = math.sqrt(2)

    # Calculate the components of the formulas
    term_M_plus_mu = M + mu
    term_M_minus_mu = M - mu
    term_div_sqrt2_1 = term_M_minus_mu / sqrt2_val
    denominator1 = term_M_plus_mu - term_div_sqrt2_1
    gamma1 = 2 / denominator1

    term_div_sqrt2_2 = term_M_minus_mu / sqrt2_val
    denominator2 = term_M_plus_mu + term_div_sqrt2_2
    gamma2 = 2 / denominator2

    # Print the step-by-step calculation for the first step size
    print("The first optimal step size (gamma_1):")
    print(f"gamma_1 = 2 / ( ({M} + {mu}) - ({M} - {mu}) / {sqrt2_val:.8f} )")
    print(f"gamma_1 = 2 / ( {term_M_plus_mu} - {term_div_sqrt2_1:.8f} )")
    print(f"gamma_1 = 2 / {denominator1:.8f}")
    print(f"gamma_1 = {gamma1:.8f}")
    
    print("")

    # Print the step-by-step calculation for the second step size
    print("The second optimal step size (gamma_2):")
    print(f"gamma_2 = 2 / ( ({M} + {mu}) + ({M} - {mu}) / {sqrt2_val:.8f} )")
    print(f"gamma_2 = 2 / ( {term_M_plus_mu} + {term_div_sqrt2_2:.8f} )")
    print(f"gamma_2 = 2 / {denominator2:.8f}")
    print(f"gamma_2 = {gamma2:.8f}")

# Execute the function
calculate_optimal_steps_for_two_step_gd()