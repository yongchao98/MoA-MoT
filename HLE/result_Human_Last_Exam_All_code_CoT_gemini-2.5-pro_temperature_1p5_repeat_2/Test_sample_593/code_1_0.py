def solve_treewidth_union():
    """
    Explains the derivation of the tight upper bound for the treewidth of the union of two graphs.
    """
    # Define symbolic variables for the explanation
    t_H = "t_H"
    t_G = "t_G"
    k = "k"

    print("This script explains the derivation of the tight upper bound for the treewidth of F.")
    print("-" * 70)

    # Step 1: Explain the upper bound derivation
    print("Step 1: Establishing an upper bound for tw(F)")
    print("Let F be the union of graphs H and G, which intersect on a set S of k vertices.")
    print("The derivation relies on the properties of clique-sums.")
    print("1. Define H' = H U K_S and G' = G U K_S, where K_S is a clique on the k shared vertices.")
    print(f"2. The treewidth of H' is bounded by tw(H') <= tw(H) + k = {t_H} + {k}.")
    print(f"3. Similarly, the treewidth of G' is bounded by tw(G') <= tw(G) + k = {t_G} + {k}.")
    print("4. F is a subgraph of H' U G', and tw(H' U G') is the maximum of tw(H') and tw(G').")
    print("This leads to the upper bound:")
    print(f"   tw(F) <= max(tw(H'), tw(G'))")
    print(f"   tw(F) <= max({t_H} + {k}, {t_G} + {k})")
    print(f"   tw(F) <= max({t_H}, {t_G}) + {k}\n")

    # Step 2: Explain the tightness proof
    print("Step 2: Proving the bound is tight")
    print("To prove tightness, we construct graphs H and G for which the equality holds.")
    print("Assume t_H >= t_G. We construct H and G such that:")
    print(f" - tw(H) = {t_H}")
    print(f" - tw(G) = {t_G}")
    print(f" - |V(H) intersect V(G)| = {k}")
    print(f" - tw(H U G) = {t_H} + {k}")
    print("The construction involves graphs where the k shared vertices (S) are an independent set connected to large cliques.")
    print("By analyzing a minor of the resulting graph F = H U G, we can show it contains a large clique.")
    print(f"This minor has a treewidth of at least {t_H} + {k}.")
    print(f"Therefore, tw(F) must be at least {t_H} + {k}.")
    print(f"Since we have tw(F) <= {t_H} + {k} (from Step 1) and tw(F) >= {t_H} + {k} (from this construction), the bound is tight.\n")

    # Step 3: State the final conclusion
    print("-" * 70)
    print("Conclusion: The tight upper bound on the treewidth of F is the following expression:")
    
    # "Output each number in the final equation!"
    # The 'numbers' here are the symbolic variables.
    print(f"max({t_H}, {t_G}) + {k}")
    print("-" * 70)

solve_treewidth_union()