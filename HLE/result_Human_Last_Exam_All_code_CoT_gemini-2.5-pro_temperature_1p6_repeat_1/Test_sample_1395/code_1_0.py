def solve_geopolitical_graph_problem():
    """
    Solves the problem by analyzing the change in the Asian country graph due to the USSR's dissolution.
    """

    # Part 1: Calculate delta_soviet, the change in the chromatic number.

    # Before the dissolution, the Asian country graph lacked a subgraph requiring 4 colors.
    # Its chromatic number, X_before, is considered to be 3.
    chi_before = 3
    print(f"Step 1: The chromatic number of the Asian subgraph before dissolution (X_before) was {chi_before}.")

    # After the dissolution, Uzbekistan and its 5 neighbors (Kazakhstan, Kyrgyzstan, Tajikistan, Afghanistan, Turkmenistan)
    # formed a W6 wheel graph, which requires 4 colors. This raised the chromatic number of the entire Asian subgraph.
    # By the Four-Color Theorem, it can be 4-colored.
    chi_after = 4
    print(f"Step 2: The chromatic number of the Asian subgraph after dissolution (X_after) became {chi_after}.")

    # The incremental change is the difference between the two.
    delta_soviet = chi_after - chi_before
    print(f"Step 3: The change in chromatic number, delta_soviet = X_after - X_before = {chi_after} - {chi_before} = {delta_soviet}.")

    # Part 2: Calculate beta, representing the change in planarity.

    # A country border graph is inherently planar. The dissolution of one large region into several
    # smaller, contiguous regions does not alter this fundamental topological property.
    # No non-planar K5 or K3,3 minors were introduced.
    print("\nStep 4: The planarity of the Asian country graph did not change after the dissolution.")

    # beta is 1 if planarity did not change, 2 if it did.
    beta = 1
    print(f"Step 5: As planarity did not change, beta is set to {beta}.")

    # Final Calculation
    final_answer = beta * delta_soviet
    print(f"\nFinal calculation: beta * delta_soviet = {beta} * {delta_soviet} = {final_answer}")


solve_geopolitical_graph_problem()