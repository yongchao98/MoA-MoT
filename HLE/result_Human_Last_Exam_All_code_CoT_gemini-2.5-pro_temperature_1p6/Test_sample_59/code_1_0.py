def print_probability_formula():
    """
    This function prints the formula for the probability of a link in a
    jointly exchangeable random graph.

    The derivation follows the standard W-graph (graphon) model, which is a
    representation for jointly exchangeable graphs.

    1.  The model assumes each node `i` has a latent variable `ξ_i` drawn
        independently from a Uniform distribution U on [0,1].
    2.  The graphon `W(u,v)` is a symmetric function that gives the probability
        of an edge between two nodes with latent variables `u` and `v`.
    3.  In the most general case, the graphon `W` itself is random, drawn
        from a distribution `F` on the space of graphons.

    The probability P(y_ij = 1) is found by averaging over all these
    random components.
    """

    # Define the components of the formula
    probability = "P(y_ij = 1)"
    expectation_operator = "E_F"  # Expectation over the random measure F
    integral_symbol = "∫"
    integration_variable_1 = "u"
    integration_variable_2 = "v"
    integration_domain = "[0,1]"
    lower_bound = 0
    upper_bound = 1
    graphon_function = f"W({integration_variable_1}, {integration_variable_2})"
    differentials = f"d{integration_variable_1} d{integration_variable_2}"

    # Construct the formula string
    # The inner part is the double integral of the graphon W(u,v) over the unit square.
    # This gives the edge density for a specific graphon W.
    inner_integral = (
        f"{integral_symbol}_{lower_bound}^{upper_bound} "
        f"{integral_symbol}_{lower_bound}^{upper_bound} "
        f"{graphon_function} {differentials}"
    )

    # The outer part is the expectation of this edge density over the distribution F of graphons.
    final_formula = f"{probability} = {expectation_operator}[ {inner_integral} ]"

    # Print the full explanation and formula
    print("For a jointly exchangeable random graph, the probability of an edge P(y_ij = 1)")
    print("is given by the expected edge density of a random graphon W.")
    print("\nHere, W is a random function drawn from a distribution F, and the expectation")
    print("of its density is taken with respect to F. The density itself is an integral")
    print("over the unit square, where latent node variables are drawn from a uniform measure U.")
    print("\nThe size of the graph, N, does not affect this probability for N >= 2.\n")
    print("The final formula is:")
    print(final_formula)


# Execute the function to print the result
print_probability_formula()