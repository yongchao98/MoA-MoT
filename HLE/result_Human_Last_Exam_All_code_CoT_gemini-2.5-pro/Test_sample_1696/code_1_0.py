import math

def solve_p_and_print_equation(n):
    """
    Calculates the closed-form formula for P(n) and prints the equation.
    The derived formula is P(n) = (2*pi)^(n*(n-1)/4) / sqrt(n!)
    
    Args:
    n (int): A positive integer.
    """
    if not isinstance(n, int) or n <= 0:
        print("Error: n must be a positive integer.")
        return

    # Calculate the components of the formula
    exponent_val = n * (n - 1) / 4
    factorial_n = math.factorial(n)
    
    # The final equation is of the form: P(n) = (2 * pi)^A / sqrt(B)
    # We output the numbers for the final equation.
    
    print(f"The closed-form formula is P(n) = (2 * pi)^(n*(n-1)/4) / sqrt(n!)")
    print(f"For n = {n}, the equation with the specific numbers is:")
    
    # Format the exponent to be an integer if it is a whole number
    if exponent_val == int(exponent_val):
        exponent_str = str(int(exponent_val))
    else:
        exponent_str = str(exponent_val)
        
    print(f"P({n}) = (2 * pi)^{exponent_str} / sqrt({factorial_n})")

# --- User can change this value ---
n_value = 5
# ---------------------------------

solve_p_and_print_equation(n_value)
