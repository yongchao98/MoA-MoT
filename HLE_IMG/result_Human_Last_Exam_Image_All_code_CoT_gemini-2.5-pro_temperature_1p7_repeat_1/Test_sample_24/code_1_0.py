def identify_manifold_from_heegaard_diagram():
    """
    Identifies a 3-manifold by analyzing the fundamental group
    derived from its Heegaard diagram.
    """

    # Step 1: The diagram is of genus 3. This gives 3 generators for the
    # fundamental group, corresponding to the beta-curves.
    num_generators = 3
    generators = [f"b{i+1}" for i in range(num_generators)]

    # Step 2: The alpha-curves provide the relators for the group presentation.
    # Each alpha-curve encircles a vertex where two beta-curves meet,
    # which implies a commutator relation.
    # The pairs of indices for the commutators are (2,3), (3,1), and (1,2).
    relator_pairs = [(2, 3), (3, 1), (1, 2)]

    print("Step 1: Determine the fundamental group presentation π₁(M).")
    print(f"The diagram has genus {num_generators}, giving {num_generators} generators: {', '.join(generators)}")
    print("\nThe relators are derived from the α-curves:")

    for i, pair in enumerate(relator_pairs):
        g1_idx, g2_idx = pair
        print(f"  - α{i+1} yields the relation: [b_{g1_idx}, b_{g2_idx}] = 1")

    # Step 3: Analyze the group defined by this presentation.
    print("\nStep 2: Simplify and identify the group.")
    print("A relation of the form [g1, g2] = 1 means g1 * g2 = g2 * g1.")
    print("The derived relations indicate that all generators commute with each other.")
    print("This defines the free abelian group on 3 generators.")
    group_name = "Z^3"
    print(f"Therefore, the fundamental group π₁(M) is isomorphic to {group_name}.")
    
    # Step 4: Identify the 3-manifold from its fundamental group.
    print(f"\nStep 3: Identify the manifold.")
    print(f"The compact, orientable 3-manifold with fundamental group {group_name} is the 3-torus.")
    
    manifold_name = "3-torus"
    manifold_symbol_base = "T"
    manifold_symbol_exponent = 3
    
    print(f"\nThe manifold is the {manifold_name}.")
    print("\nFinal symbolic equation for the manifold:")
    print(f"M = {manifold_symbol_base}^{manifold_symbol_exponent}")
    print("\nThe numbers in the symbolic equation are:")
    print(f"  - Exponent: {manifold_symbol_exponent}")

identify_manifold_from_heegaard_diagram()