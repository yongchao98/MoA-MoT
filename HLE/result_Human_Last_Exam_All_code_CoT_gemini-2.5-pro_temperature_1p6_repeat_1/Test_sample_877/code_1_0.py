import math

def print_final_equation():
    """
    This function prints the derived equation for h(b(0)) by constructing it
    from its numerical components, as requested.
    """
    coeff_b_sq = 4
    coeff_b = -6
    const = 2
    coeff_bln = 2
    factor_ln = 2
    
    # The variable in the function is b(0)
    var = 'b(0)'
    
    print(f"The function h({var}) is determined to be:")
    
    # Print the equation piece by piece to show each number
    # The str.format method is used to control spacing and signs.
    # We build the string expression for the equation.
    
    equation_str = "{0}*{1}^2 ".format(coeff_b_sq, var)
    
    if coeff_b < 0:
        equation_str += "- {0}*{1} ".format(-coeff_b, var)
    else:
        equation_str += "+ {0}*{1} ".format(coeff_b, var)

    if const < 0:
        equation_str += "- {0} ".format(-const)
    else:
        equation_str += "+ {0} ".format(const)

    if coeff_bln < 0:
        equation_str += "- {0}*{1}*ln({2}*{1})".format(-coeff_bln, var, factor_ln)
    else:
        equation_str += "+ {0}*{1}*ln({2}*{1})".format(coeff_bln, var, factor_ln)
        
    print("h({0}) = {1}".format(var, equation_str))

# Execute the function to print the final equation
print_final_equation()
