def print_solution():
    """
    This function prints the general solution of the given differential equation.
    The solution was found through the derivation explained above.
    The equation is y**2 = C*x**2 + C**2 + C + 9, where C is an arbitrary constant.
    To satisfy the prompt's request, we will format the output to clearly show the numbers
    in one of its standard forms, like: 1*y**2 - C*x**2 - 1*C**2 - 1*C - 9 = 0.
    """
    
    # Let C represent the arbitrary constant.
    C = 'C'
    
    # The general solution equation can be written as:
    # 1*y**2 - C*x**2 - C**2 - C - 9 = 0
    
    # Coefficients and terms
    coeff_y_sq = 1
    term_x_sq = f"-{C}*x**2"
    term_C_sq = f"-{1}*{C}**2"
    term_C = f"-{1}*{C}"
    constant_term = -9
    
    print("The general solution of the differential equation is:")
    print(f"({coeff_y_sq})*y**2 + ({term_x_sq.split('*')[0]})*x**2 + ({term_C_sq.split('*')[0]})*{C}**2 + ({term_C.split('*')[0]})*{C} + ({constant_term}) = 0")
    print("where C is an arbitrary constant.")
    
    print("\nA more common way to write this is:")
    print(f"y**2 = {C}*x**2 + {C}**2 + {C} + {9}")

print_solution()