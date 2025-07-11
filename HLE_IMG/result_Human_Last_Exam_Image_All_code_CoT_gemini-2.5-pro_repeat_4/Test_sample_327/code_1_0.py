def solve_graphene_puzzle():
    """
    Analyzes the band structures of three graphene nanoribbons and generates a combined classification string.
    """

    # --- Analysis of Ribbon 1 ---
    # Edge: Armchair (A) - Symmetric band structure, gap at k=0.
    # Width: N=8 - 8 bands counted above the Fermi level (E=0).
    # Band: Metallic (0) - Conduction and valence bands touch at E=0.
    edge1 = "A"
    width1 = 8
    band1 = 0
    classification1 = f"{edge1}{width1}{band1}"

    # --- Analysis of Ribbon 2 ---
    # Edge: Zigzag (Z) - Characteristic flat bands near the Fermi level.
    # Width: N=7 - 7 bands counted above the Fermi level at k=0.
    # Band: Metallic (0) - Bands cross the Fermi level.
    edge2 = "Z"
    width2 = 7
    band2 = 0
    classification2 = f"{edge2}{width2}{band2}"

    # --- Analysis of Ribbon 3 ---
    # Edge: Armchair (A) - Symmetric band structure, gap at k=0.
    # Width: N=7 - 7 bands counted above the Fermi level.
    # Band: Semiconducting (1) - Clear energy gap around the Fermi level.
    edge3 = "A"
    width3 = 7
    band3 = 1
    classification3 = f"{edge3}{width3}{band3}"

    # Concatenate all classifications into the final result string
    final_result = classification1 + classification2 + classification3

    print(final_result)

solve_graphene_puzzle()