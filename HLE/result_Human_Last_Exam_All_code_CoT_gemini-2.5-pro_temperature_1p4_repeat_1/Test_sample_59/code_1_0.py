def solve_probability_formula():
    """
    This script symbolically constructs and prints the formula for the probability
    of an edge in a jointly exchangeable random graph.
    """
    
    # --- Define the symbols used in the formula ---
    
    # The probability we want to calculate
    probability_term = "P(y_ij = 1)"
    
    # The expectation over the random measure F, from which the graphon W is drawn
    expectation_F = "E_F"
    
    # The integral symbol and its limits
    integral_symbol = "∫"
    lower_limit = 0
    upper_limit = 1
    
    # The graphon function W(u,v) and the differentials
    graphon_function = "W(u,v)"
    differentials = "du dv"

    # --- Construct and print the final equation ---
    
    # The inner part: the double integral of the graphon
    integral_part = (f"{integral_symbol}_{lower_limit}^{upper_limit} "
                     f"{integral_symbol}_{lower_limit}^{upper_limit} "
                     f"{graphon_function} {differentials}")

    # The full formula, including the expectation over F
    final_equation = f"{probability_term} = {expectation_F}[ {integral_part} ]"

    print("The formula for the probability of a link y_ij in a jointly exchangeable random graph is:")
    print(final_equation)
    
    print("\n--- Components of the Equation ---")
    print(f"The value representing an existing link (y_ij) is: 1")
    print(f"The lower limit for the integrals is: {lower_limit}")
    print(f"The upper limit for the integrals is: {upper_limit}")
    print("\nWhere:")
    print("P(y_ij = 1) is the probability of an edge between nodes i and j.")
    print("E_F is the expectation over the distribution F of random graphons.")
    print("W(u,v) is the graphon function, which gives the conditional edge probability.")
    print("u and v are latent variables for nodes i and j, integrated over the interval [0,1].")

solve_probability_formula()