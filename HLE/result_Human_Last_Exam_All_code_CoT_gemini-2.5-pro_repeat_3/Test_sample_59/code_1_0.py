import sympy

def solve_graph_probability():
    """
    Explains and derives the formula for the probability of a link in a
    jointly exchangeable random graph.
    """

    print("Deriving the probability of a link y_ij in a jointly exchangeable random graph.")
    print("-" * 70)

    # Step 1: Explain the underlying model (Graphon model)
    print("\nStep 1: The Model")
    print("According to the Aldous-Hoover theorem for exchangeable arrays, a jointly exchangeable graph can be generated from a model called a 'graphon'.")
    print("1. For each node `i` in the graph, we assign a latent variable `u_i` drawn independently from a Uniform distribution on [0,1].")
    print("   u_i ~ U[0,1]")
    print("2. A symmetric, measurable function W(x, y), the 'graphon', maps the unit square [0,1]x[0,1] to a probability [0,1].")
    print("3. The probability of a link `y_ij` between nodes `i` and `j` is determined by their latent variables through the graphon W.")
    print("\nP(y_ij = 1 | u_i, u_j) = W(u_i, u_j)")


    # Step 2: Derive the marginal probability
    print("\nStep 2: Calculating the Marginal Probability")
    print("The overall probability of a link, P(y_ij = 1), is the expectation of the conditional probability over the distribution of the latent variables u_i and u_j.")
    print("\nP(y_ij = 1) = E[ W(u_i, u_j) ]")
    print("\nSince u_i and u_j are independent and identically distributed as U[0,1], this expectation is calculated by integrating over the unit square.")

    # Step 3: Present the final equation
    # Using symbolic representation for clarity
    x, y = sympy.symbols('x y')
    W = sympy.Function('W')
    
    lower_bound_x = 0
    upper_bound_x = 1
    lower_bound_y = 0
    upper_bound_y = 1
    integrand = f"{W(x, y)}"
    differentials = "dx dy"

    print("\n" + "="*30 + " FINAL EQUATION " + "="*30)
    # We print each number and symbol in the final equation as requested.
    print(f"P(y_ij = 1) = integral from {lower_bound_y} to {upper_bound_y} ( integral from {lower_bound_x} to {upper_bound_x} of {integrand} ) {differentials}")
    print("-" * 70)
    print("\nThis probability is a constant, often denoted by 'p', which represents the expected edge density of the graph.")
    print("It depends on the specific choice of the graphon function W(x, y). Without a specific W, the probability is represented by this integral.")

solve_graph_probability()