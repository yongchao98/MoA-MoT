def count_components():
    """
    This function analyzes the number of connected components of the space X
    after the origin is removed.
    """

    # The original space X is connected because all line segments meet at the origin.
    # When the origin is removed, the single connection point is gone.

    # The space splits into multiple disjoint connected components.

    # Component 1: The segment L from (1,0) to the origin, with the origin removed.
    base_component = "The segment from (1, 0) to (0,0), with (0,0) removed."
    
    # The other components are formed from the segments L_n.
    # There is one such component for each positive integer n.
    # L_n is the segment from (1, 1/n) to the origin.

    print("The space is broken into the following connected components:")
    print("1. {}".format(base_component))
    print("2. An infinite series of components, one for each positive integer n, where the nth component is:")
    print("   The segment from (1, 1/n) to (0,0), with (0,0) removed.\n")
    
    print("Let's list the first few of these components to see the pattern:")
    
    components = []
    components.append({
        "index": "base", 
        "description": "from p=(1,0)"
    })
    
    # We list components for n = 1, 2, 3, ...
    for n in range(1, 6):
        components.append({
            "index": n,
            "description": "from p_{}=({}, 1/{})".format(n, 1, n)
        })

    for comp in components:
        print("Component related to p_{}: {}".format(comp["index"], comp["description"]))

    print("...")
    print("\nThis sequence of components continues for all integers n >= 1.")
    print("Thus, we have 1 component (from L) + an infinite number of components (from each L_n).")
    print("The total number of connected components is infinite.")

if __name__ == "__main__":
    count_components()
