def print_probability_formula():
    """
    This function prints the formula for the probability of drawing a link
    in a jointly exchangeable random graph.
    """
    
    # Define the components of the equation as strings
    probability_term = "P(y_ij = 1)"
    equals_sign = "="
    integral_symbol = "\u222B" # Unicode for the integral sign
    
    # The integral part represents averaging over the latent variables u and v,
    # which are drawn from a uniform distribution on [0,1].
    integral_expression = f"{integral_symbol}[0,1] {integral_symbol}[0,1] f(u,v) du dv"
    
    # The full equation
    # P(y_ij = 1) = integral from 0 to 1, integral from 0 to 1 of f(u,v) du dv
    final_equation = f"{probability_term} {equals_sign} {integral_expression}"

    print("The probability P(y_ij = 1) to draw a link between nodes i and j is given by the integral of the graphon function f(u,v) over the unit square.")
    print("This represents the expected density of edges in the graph.")
    print("\nThe final formula is:")
    print(final_equation)

print_probability_formula()