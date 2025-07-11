import sys

def solve():
    """
    This function explains the reasoning and prints the final answer.
    """

    # Step 1: Establish a lower bound for the number of components.
    # A fundamental theorem in hyperspace theory (G. Beer, 1993) states that for any metric space X,
    # the hyperspace CL(X) with the Wijsman topology is connected if and only if X is connected.
    # The space X in the problem is stated to be totally-disconnected and have infinitely many points.
    # A totally-disconnected space with more than one point is, by definition, not connected.
    # Therefore, CL(X) cannot be connected, which means its number of connected components must be at least 2.

    # Step 2: Show that this lower bound is achievable.
    # We need to find an X that satisfies the problem's conditions and results in CL(X) having exactly 2 components.
    # A theorem by G. Vidossich (1970) provides the necessary tool. It states that for a non-compact,
    # locally compact, Polish (complete and separable) metric space X, the hyperspace CL(X)
    # has exactly two connected components. These components are the set of all non-empty compact
    # subsets, K(X), and the set of all non-empty non-compact closed subsets, NCL(X).

    # Step 3: Identify a space X with all the required properties.
    # We need an X that is:
    # 1. An infinite, totally-disconnected ultrametric space (from the problem statement).
    # 2. A non-compact, locally compact, Polish space (to apply Vidossich's theorem).
    # The field of p-adic numbers, Q_p (for any prime p), is a perfect example.
    # - It is an ultrametric space.
    # - It is totally disconnected and infinite.
    # - It is non-compact (it is unbounded) but locally compact (closed balls are compact).
    # - It is a Polish space (it's the completion of the rationals Q under the p-adic metric).
    # Since Q_p meets all these criteria, we can apply Vidossich's theorem.
    # Thus, for X = Q_p, CL(X) has exactly 2 connected components.

    # Step 4: Conclusion.
    # We've shown that the number of components must be at least 2, and that it can be exactly 2.
    # Therefore, the smallest possible number of connected components is 2.

    smallest_number_of_components = 2

    print("The final equation is: Smallest possible number = 2")
    # To follow the instruction "output each number in the final equation!", we print the number itself.
    print(smallest_number_of_components)

solve()
