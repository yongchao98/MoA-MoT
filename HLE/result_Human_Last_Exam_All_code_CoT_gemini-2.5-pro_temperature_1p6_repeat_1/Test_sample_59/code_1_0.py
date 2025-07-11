import sympy

def solve_probability_formula():
    """
    This function symbolically computes and prints the formula for the probability 
    of an edge in a jointly exchangeable random graph.
    """

    # 1. Define the symbolic variables for the integral and the functions.
    # u and v are the integration variables, representing the latent positions of nodes.
    u, v = sympy.symbols('u v')
    
    # W is the graphon function W(u, v), which gives the probability of an edge 
    # for given latent positions u and v.
    W = sympy.Function('W')
    
    # E represents the expectation operator. It is needed because the graphon W
    # itself is random, drawn from a distribution determined by the random measure F.
    E = sympy.Function('E')

    # 2. Construct the formula.
    # The probability of an edge for a *specific* graphon W is the integral of W(u,v)
    # over the unit square [0,1] x [0,1]. The limits of integration are 0 and 1.
    integral_of_W = sympy.Integral(W(u, v), (v, 0, 1), (u, 0, 1))
    
    # The final probability P(y_ij = 1) is the expectation of this integral over
    # the distribution of all possible graphons W.
    final_probability_formula = E(integral_of_W)

    # 3. Print the result in a readable format.
    print("The general formula for the probability of a link y_{ij} in a jointly exchangeable random graph is:")
    
    # sympy.pretty_print() renders the formula in a more mathematical style.
    sympy.pretty_print(final_probability_formula)

if __name__ == "__main__":
    solve_probability_formula()