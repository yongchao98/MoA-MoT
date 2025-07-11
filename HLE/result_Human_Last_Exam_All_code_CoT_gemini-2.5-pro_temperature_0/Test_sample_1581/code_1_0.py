def solve():
    """
    This function solves the topological problem posed by the user.

    The problem asks for the number of distinct homeomorphism classes for a compact
    connected metric space X, given that its n-th configuration space F_n(X) is
    disconnected for some n >= 2.

    The reasoning is as follows:
    1.  If X is an arc (homeomorphic to [0,1]), its total ordering causes F_n(X)
        to be disconnected for all n >= 2. The space of configurations splits
        into n! components based on the ordering of the points, and no path
        can exist between these components. So, the class of arcs is a solution.

    2.  Conversely, a key theorem in topology states that for a broad class of
        spaces (including all Peano continua), F_n(X) is connected if and only
        if X is NOT an arc. This result also extends to common non-locally
        connected spaces like the topologist's sine curve.

    3.  Therefore, the condition that F_n(X) is disconnected essentially forces X
        to be an arc.

    4.  All arcs are homeomorphic to each other. They form a single
        homeomorphism class.

    5.  Thus, there is only one such homeomorphism class.
    """
    
    # The number of distinct homeomorphism classes.
    number_of_classes = 1
    
    # The problem asks to output the numbers in the final equation.
    # Here, the final answer is simply the number itself.
    print(f"The number of distinct homeomorphism classes is: {number_of_classes}")

solve()