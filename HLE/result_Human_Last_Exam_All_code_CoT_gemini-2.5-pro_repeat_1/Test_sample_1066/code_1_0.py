import math

def print_lower_bound_equation():
    """
    This function prints the derived lower bound for the expected detection score E[S].
    The formula is E[S] >= n * alpha * (pi^2/6 - 1) + n.
    The code will calculate the numerical value of the constant factor.
    """
    
    # Calculate the constant part of the expression
    pi_squared_over_6 = math.pi**2 / 6
    constant_factor = pi_squared_over_6 - 1
    
    # The final equation string
    equation = f"E[S] >= n * alpha * ({pi_squared_over_6:.4f} - 1) + n"
    final_equation = f"E[S] >= n * alpha * {constant_factor:.4f} + n"
    
    print("The derived lower bound for the expected score E[S] is:")
    print(equation)
    print("Which simplifies to:")
    print(final_equation)

print_lower_bound_equation()