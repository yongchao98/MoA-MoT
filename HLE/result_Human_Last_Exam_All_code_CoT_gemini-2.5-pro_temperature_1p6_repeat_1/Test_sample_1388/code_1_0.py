def get_H_t_expression():
    """
    This function prints the explicit form of H(t) based on the derivation.
    The final equation includes the numbers 1 and 2 explicitly as requested.
    """
    
    # The numerator and denominator in the exponent.
    numerator = 1
    denominator = 2
    
    # Definition of h(t) as an integral expression.
    h_t_definition = "Integral_{s=0 to t} ||∂x u(·, s)||_{L∞} ds"
    
    # Construct the final expression for H(t).
    H_t_expression = f"exp( ({numerator}/{denominator}) * ( {h_t_definition} ) )"
    
    print("The explicit form of H(t) is:")
    print(H_t_expression)

get_H_t_expression()