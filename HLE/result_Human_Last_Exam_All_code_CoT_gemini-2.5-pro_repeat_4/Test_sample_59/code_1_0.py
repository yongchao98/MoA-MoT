def explain_link_probability():
    """
    Explains and derives the probability of a link in a jointly exchangeable random graph.
    """
    print("Derivation of the Probability to Draw a Link in a Jointly Exchangeable Graph")
    print("==============================================================================")
    print("\nA jointly exchangeable random graph has a structure whose probability is invariant")
    print("to the relabeling of its N nodes. The probability of a link y_ij between any")
    print("two nodes i and j is determined by the Aldous-Hoover representation theorem.")
    print("\nStep 1: Latent Variables")
    print("--------------------------")
    print("Each node 'i' is assigned a latent (hidden) variable, alpha_i, drawn independently")
    print("from a Uniform distribution on the interval [0, 1].")
    print("  alpha_i ~ U(0, 1)")

    print("\nStep 2: The Graphon Function f(u, v)")
    print("-----------------------------------")
    print("A symmetric function f(u, v), called a 'graphon', defines the probability of a")
    print("link between two nodes based on their latent variables.")

    print("\nStep 3: Conditional Link Probability")
    print("------------------------------------")
    print("The probability of a link y_ij = 1, *given* the latent variables alpha_i and alpha_j, is:")
    print("  P(y_ij = 1 | alpha_i, alpha_j) = f(alpha_i, alpha_j)")

    print("\nStep 4: Unconditional Link Probability")
    print("--------------------------------------")
    print("To find the unconditional probability P(y_ij = 1), we average the conditional")
    print("probability over all possible values of alpha_i and alpha_j. Since they are drawn")
    print("from a Uniform distribution, this average is a double integral over the unit square [0,1]x[0,1].")

    print("\nFinal Equation:")
    print("P(y_ij = 1) = integral from 0 to 1 of ( integral from 0 to 1 of f(u, v) du ) dv")
    print("\nThis integral represents the expected edge density of the graph. Due to the exchangeability,")
    print("the probability of a link is the same for any pair of distinct nodes (i, j).")

explain_link_probability()

# The final answer is the formula derived above.
# It's not a single number because the function f is not specified.
# The formula itself is the answer.
final_answer = "P(y_ij = 1) = integral from 0 to 1 ( integral from 0 to 1 ( f(u, v) du ) dv )"