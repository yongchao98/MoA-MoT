import sys

def solve_probability():
    """
    Derives and explains the probability of drawing a link in a jointly
    exchangeable random graph.
    """
    print("Derivation of the Link Probability in a Jointly Exchangeable Random Graph")
    print("-" * 70)

    # Step 1: Explain Joint Exchangeability
    print("1. Understanding Joint Exchangeability")
    print("A random graph, defined by its adjacency matrix Y = (y_ij), is jointly exchangeable")
    print("if its probability distribution is unchanged by any relabeling (permutation) of its nodes.")
    print("A direct consequence of this symmetry is that the probability of an edge existing")
    print("between any two distinct nodes i and j is the same for all such pairs.")
    print("\n  P(y_ij = 1) = p   (for all i != j)")
    print("\nOur goal is to find an expression for this constant probability, p.")
    print("\n")

    # Step 2: The Aldous-Hoover Representation Theorem
    print("2. The Aldous-Hoover Representation Theorem")
    print("This fundamental theorem provides a way to construct any jointly exchangeable random graph.")
    print("It states that such a graph can be generated as follows:")
    print("  - To each node i, we assign a latent variable xi_i, drawn independently from a")
    print("    Uniform(0,1) distribution. This corresponds to the 'U and uniform measure' in the prompt.")
    print("  - There exists a symmetric, measurable function W: [0,1] x [0,1] -> [0,1], called a graphon.")
    print("    This graphon W can be random, generated from a process described by 'F a random measure'.")
    print("  - The probability of an edge between nodes i and j, given their latent variables and W, is:")
    print("\n    P(y_ij = 1 | xi_i, xi_j, W) = W(xi_i, xi_j)")
    print("\n")

    # Step 3: Deriving the Unconditional Probability
    print("3. Calculation of the Probability")
    print("To find the unconditional probability p = P(y_ij = 1), we use the law of total probability.")
    print("We must average over all possible values of the latent variables xi_i and xi_j, and also")
    print("over the randomness in the graphon W itself.")
    print("\n  p = E[ W(xi_i, xi_j) ]")
    print("\nSince xi_i and xi_j are independent and uniformly distributed on [0,1], we can express this")
    print("expectation as an integral. The outer expectation E_W averages over the choice of graphon W.")
    print("\n  p = E_W [ integral from 0 to 1 of integral from 0 to 1 of W(u,v) du dv ]")
    print("\nThe double integral inside the expectation, ∫∫W(u,v)dudv, represents the edge density of a")
    print("specific graphon W. Therefore, the link probability is the expected edge density of the graph.")
    print("\n")

    # Step 4: Printing the Final Equation
    print("4. The Final Equation for the Link Probability")
    print("The general expression for the probability to draw a link y_ij is:")
    print("\n  P(y_ij = 1) = E[∫_0^1 ∫_0^1 W(u,v) du dv]")
    print("\nPrinting each symbol of the final equation as requested:")
    # Using sys.stdout.write to avoid newlines and spaces for more control
    sys.stdout.write("P ( y_ij = 1 ) = E [ ∫_0^1 ∫_0^1 W ( u , v ) du dv ]\n")


solve_probability()