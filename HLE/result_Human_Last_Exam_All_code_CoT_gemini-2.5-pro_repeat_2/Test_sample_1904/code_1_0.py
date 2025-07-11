def solve():
    """
    This function solves the problem about the number of connected components of CL(X).
    """
    
    # Step 1: A known theorem states that CL(X) is connected if and only if X is connected.
    # X is totally-disconnected and has infinitely many points, so it is not connected.
    # Therefore, CL(X) is not connected, and the number of components must be at least 2.
    lower_bound = 2
    
    # Step 2: Researchers have shown that for a non-compact, unbounded ultrametric space X,
    # the connected components of CL(X) are determined by the 'ends' of the space.
    # One component consists of all metrically bounded closed sets.
    # The other components consist of classes of unbounded sets, corresponding to the ends of X.
    
    # Step 3: It is possible to construct an infinite, totally-disconnected ultrametric space X
    # that has only a single 'end'.
    
    # Step 4: For such a space with a single end, there is one component for bounded sets
    # and one component for all unbounded sets.
    # This leads to a total of 1 + 1 = 2 components.
    
    # Step 5: Since the number of components is at least 2, and we can achieve exactly 2,
    # the smallest possible number is 2.
    smallest_possible_number = 2
    
    print(f"The smallest possible number of connected components of CL(X) is {smallest_possible_number}.")

solve()