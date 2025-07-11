def explain_exchangeable_graph_probability():
    """
    This function explains and prints the formula for the probability of an edge
    in a jointly exchangeable random graph.
    """

    print("The probability P(y_ij = 1) for a link in a jointly exchangeable random graph is derived from the Aldous-Hoover representation theorem for such graphs.")
    print("This theorem states that the graph can be generated from an underlying object called a graphon, W(u,v).")
    print("-" * 70)

    # Symbolically define the components of the final formula
    final_equation_parts = {
        "The probability of an edge": "P(y_ij = 1)",
        "is equal to": "=",
        "The expectation over the random graphon W": "E_W",
        "of the following quantity": "[",
        "The integral over latent variable u from 0 to 1": "Integral(u=0 to 1)",
        "of the integral over latent variable v from 0 to 1": "Integral(v=0 to 1)",
        "of the graphon function W(u,v)": "W(u,v)",
        "with respect to u and v": "du dv",
        "end of expectation": "]"
    }

    print("The final formula is constructed by averaging over all sources of randomness:")
    print("P(y_ij = 1) = E_W [ Integral(from u=0 to 1) Integral(from v=0 to 1) W(u,v) du dv ]")

    print("\n--- Components of the Final Equation ---")
    # We output each component of the formula as requested.
    for description, symbol in final_equation_parts.items():
        print(f"'{symbol}'  ({description})")

    print("-" * 70)
    print("\nExplanation of components:")
    print("  - W(u,v): The random graphon, a function representing the edge probability between nodes with latent positions u and v.")
    print("  - E_W[...]: The expectation over the random choice of the graphon W (drawn from what the problem calls a random measure 'F').")
    print("  - Integral(...): A double integral that averages over all possible latent positions u and v, which are drawn from a Uniform measure 'U' on [0,1].")

explain_exchangeable_graph_probability()