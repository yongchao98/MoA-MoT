import sympy

def solve_probability_exchangeable_graph():
    """
    Derives the probability of an edge in a jointly exchangeable random graph.

    The derivation is presented symbolically using the sympy library.
    """

    # Define the mathematical symbols and functions required for the derivation.
    # u and v represent the latent variables for nodes i and j. They are drawn from a
    # Uniform measure U on the interval [0,1].
    u, v = sympy.symbols('u v')

    # F represents the random measure that parameterizes the graph's structure.
    # W_F(u, v) is the graphon function. It's a random function because it depends on F.
    # It gives the probability of an edge between nodes with latent variables u and v.
    W_F = sympy.Function('W')(u, v)

    # --- Step-by-step derivation ---

    print("Derivation of the probability P(y_ij = 1):\n")

    print("Step 1: Conditional Probability given latent variables and the random measure.")
    print("According to the Aldous-Hoover theorem for exchangeable graphs (graphon model),")
    print("the probability of an edge y_ij between nodes i and j, conditional on their")
    print(f"latent variables U_i=u, U_j=v and the random measure F, is given by the graphon function W(u,v):")
    # Using a generic W for simplicity in the printout
    print(f"P(y_ij = 1 | U_i=u, U_j=v, F) = {sympy.Function('W')(u,v)}\n")

    print("Step 2: Averaging over the latent variables.")
    print("To find the probability for a fixed random measure F, we average over all possible values of")
    print("U_i and U_j. Since they are drawn independently from the Uniform distribution on [0,1],")
    print("we integrate W(u,v) over the unit square.")
    
    # This integral represents the expected edge density for a specific realization of the graphon W.
    edge_density_for_fixed_F = sympy.Integral(W_F, (u, 0, 1), (v, 0, 1))
    
    print("The resulting probability, conditional only on F, is rho(F):")
    print(f"rho(F) = E[y_ij | F] = {edge_density_for_fixed_F}\n")

    print("Step 3: Averaging over the random measure F.")
    print("The final, unconditional probability is the expectation of rho(F) over the")
    print("distribution of F. Let E_F[...] denote this expectation.")
    
    # We use sympy's expectation symbol E to represent the final step.
    # For a cleaner final expression, we represent the integral as rho(F).
    rho_F = sympy.Function('rho')(sympy.Symbol('F'))
    final_probability_expression = sympy.E(rho_F)
    
    print(f"P(y_ij = 1) = E_F[rho(F)] = {final_probability_expression}\n")

    print("--- Final Equation ---")
    print("Substituting the definition of rho(F) back into the expression gives the full formula.")
    print("Each number and symbol in the final equation is shown below:")
    
    # Final equation with all components visible
    final_equation_str = f"P(y_ij = 1) = E_F[ Integral({W_F}, (u, {0}, {1}), (v, {0}, {1})) ]"
    print(final_equation_str)
    
    print("\nIn words: The probability to draw a link is the expected edge density, where the")
    print("expectation is taken over the distribution of random graphons.")

solve_probability_exchangeable_graph()
<<<E_F[ Integral(W(u,v), (u, 0, 1), (v, 0, 1)) ]>>>