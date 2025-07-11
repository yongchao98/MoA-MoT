import numpy as np

def print_lower_bound_equation():
    """
    This function prints the derived lower bound for E[S].
    """
    
    # Define the variables for the equation string
    E_S = "E[S]"
    n = "n"
    alpha = "alpha"
    pi = "pi"
    
    # The lower bound E[S] >= n*alpha + n*(1 - pi^2/6)
    bound_formula = f"{E_S} >= {n}*{alpha} + {n}*(1 - {pi}**2 / 6)"
    
    print("The derived lower bound on the expected score E[S] is:")
    print(bound_formula)
    
    print("\nWhere the numerical values for the constants are approximately:")
    pi_val = np.pi
    constant_val = 1 - (pi_val**2) / 6
    print(f"pi = {pi_val:.6f}")
    print(f"1 - pi**2 / 6 = {constant_val:.6f}")

    # This part prints the final equation with numbers as requested by the prompt.
    print("\nFinal equation with each number printed separately:")
    n_val = "n"
    alpha_val = "alpha"
    one_val = 1
    six_val = 6
    pi_symbol = "pi"

    print(f"E[S] >= {n_val} * {alpha_val} + {n_val} * ({one_val} - {pi_symbol}**2 / {six_val})")


print_lower_bound_equation()
