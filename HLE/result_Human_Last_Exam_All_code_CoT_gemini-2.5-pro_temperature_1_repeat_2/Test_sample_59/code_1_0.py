def explain_exchangeable_graph_probability():
    """
    This function explains and prints the formula for the probability of an edge
    in a jointly exchangeable random graph.
    """

    # Define the components of the formula as symbolic strings.
    prob_of_edge = "P(y_ij = 1)"
    expectation_operator = "E"
    integral_lower_limit = "0"
    integral_upper_limit = "1"
    integral_symbol = "∫"
    graphon_function = "W(u, v)"
    differentials = "du dv"

    # Assemble the final formula string.
    # The formula represents the expectation (over the random choice of graphon W)
    # of the integral of the graphon over the unit square.
    final_formula = (
        f"{prob_of_edge} = {expectation_operator}["
        f"{integral_symbol}_{integral_lower_limit}^{integral_upper_limit} "
        f"{integral_symbol}_{integral_lower_limit}^{integral_upper_limit} "
        f"{graphon_function} {differentials}]"
    )

    print("For a jointly exchangeable random graph, the probability of an edge y_ij is given by the Aldous-Hoover theorem.")
    print("The formula for this probability is the expected edge density:\n")

    # Print the final formula, showing each component including the numbers for the limits.
    print(final_formula)

    print("\nWhere:")
    print(f" - {prob_of_edge}: The probability of an edge between nodes i and j.")
    print(f" - {expectation_operator}: The expectation over the random measure F, which defines the random graphon W.")
    print(f" - {integral_symbol}_{integral_lower_limit}^{integral_upper_limit}: The definite integral from {integral_lower_limit} to {integral_upper_limit}.")
    print(f" - {graphon_function}: The graphon, a symmetric function defining edge probabilities from latent variables u and v.")
    print(f" - u, v: Latent variables drawn independently from a Uniform({integral_lower_limit}, {integral_upper_limit}) distribution.")

# Execute the function to print the explanation and formula.
explain_exchangeable_graph_probability()