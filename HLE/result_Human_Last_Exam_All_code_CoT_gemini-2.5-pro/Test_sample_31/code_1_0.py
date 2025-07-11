def calculate_min_pll_stickers():
    """
    Calculates the minimum number of non-top-facing stickers needed
    to identify a PLL case on a 3x3 Rubik's Cube.
    """

    # --- Step 1: Determine stickers needed for the corner permutation ---
    
    # There are 4 corners in the last layer.
    num_corners = 4
    # To determine the permutation of N items, we only need to locate N-1 of them.
    # The position of the last item is then fixed.
    corners_to_identify = num_corners - 1
    
    # Each corner piece has 2 non-top-facing (side) stickers.
    stickers_per_corner = 2
    
    # Calculate the number of stickers needed to identify the corner permutation.
    corner_stickers_needed = corners_to_identify * stickers_per_corner
    
    print(f"To determine the corner permutation of {num_corners} pieces, we need to identify {corners_to_identify}.")
    print(f"Stickers needed for corners = {corners_to_identify} corners * {stickers_per_corner} stickers/corner = {corner_stickers_needed}")
    print("-" * 20)

    # --- Step 2: Determine stickers needed for the edge permutation using parity ---

    # There are 4 edges in the last layer.
    num_edges = 4
    # By determining the full corner permutation, we know its parity.
    # The edge permutation must have the same parity.
    # We only need to identify N-2 edges. This leaves the last 2 edges in the last 2 slots.
    # This creates two possibilities for the edge permutation (swapped or not swapped),
    # which have opposite parities. The parity rule resolves this ambiguity.
    edges_to_identify = num_edges - 2
    
    # Each edge piece has 1 non-top-facing (side) sticker.
    stickers_per_edge = 1
    
    # Calculate the number of stickers needed to constrain the edge permutation.
    edge_stickers_needed = edges_to_identify * stickers_per_edge
    
    print("The cube's parity rule states: parity(corners) must equal parity(edges).")
    print(f"Knowing the corner permutation tells us the required parity for the edges.")
    print(f"We only need to identify {edges_to_identify} edges to leave a single solvable ambiguity.")
    print(f"Stickers needed for edges = {edges_to_identify} edges * {stickers_per_edge} sticker/edge = {edge_stickers_needed}")
    print("-" * 20)
    
    # --- Step 3: Calculate the total ---
    
    total_stickers = corner_stickers_needed + edge_stickers_needed
    
    print("Total minimum stickers required = (stickers for corners) + (stickers for edges)")
    print(f"Total minimum stickers required = {corner_stickers_needed} + {edge_stickers_needed} = {total_stickers}")
    
    return total_stickers

if __name__ == '__main__':
    calculate_min_pll_stickers()