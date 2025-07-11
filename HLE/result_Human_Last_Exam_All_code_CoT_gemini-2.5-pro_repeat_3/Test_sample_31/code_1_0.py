import sys

def solve_pll_stickers():
    """
    Calculates the minimum number of non-top-facing stickers
    needed to be seen to identify any PLL case.
    """
    
    # --- Step 1: Determine the corner permutation ---
    # To determine the permutation of 4 unique pieces, we only need to
    # find the identity and position of 3 of them. The 4th is then known by elimination.
    corners_to_identify = 3
    
    # To uniquely identify a corner piece, we must see its 2 side stickers
    # (e.g., the Blue and Red stickers on the Blue-Red-Yellow corner).
    stickers_per_corner = 2
    
    # Calculate the number of stickers needed to determine the full corner permutation.
    corner_stickers_needed = corners_to_identify * stickers_per_corner
    
    # --- Step 2: Determine the edge permutation using the parity rule ---
    # An important rule is that the parity of the corner permutation (even/odd)
    # must match the parity of the edge permutation. Once we know the corner
    # permutation, we know the required parity for the edges.
    
    # We only need to identify 2 of the 4 edges. The remaining two can either be
    # in their correct places or swapped. These two possibilities have different parities.
    # The parity rule allows us to determine which case is correct.
    edges_to_identify = 2
    
    # To uniquely identify an edge piece, we must see its 1 side sticker.
    stickers_per_edge = 1
    
    # Calculate stickers needed to determine the edge permutation, given the corner info.
    edge_stickers_needed = edges_to_identify * stickers_per_edge
    
    # --- Step 3: Calculate the total ---
    total_stickers = corner_stickers_needed + edge_stickers_needed

    print("To fully identify a PLL case, we must determine the permutation of the 4 corner and 4 edge pieces.")
    print("Here is the breakdown of the minimum number of stickers that must be seen:")
    print("\n1. Corner Piece Identification:")
    print(f"   - To know the permutation of 4 corners, we only need to locate {corners_to_identify}.")
    print(f"   - Each corner requires seeing {stickers_per_corner} stickers to be identified.")
    print(f"   - Stickers for corners: {corners_to_identify} * {stickers_per_corner} = {corner_stickers_needed}")
    
    print("\n2. Edge Piece Identification (using the Parity Rule):")
    print(f"   - After finding the corner permutation, we know the required parity for the edge permutation.")
    print(f"   - We only need to locate {edges_to_identify} edges. The locations of the other two are determined by the parity rule.")
    print(f"   - Each edge requires seeing {stickers_per_edge} sticker to be identified.")
    print(f"   - Stickers for edges: {edges_to_identify} * {stickers_per_edge} = {edge_stickers_needed}")
          
    print("\n3. Final Calculation:")
    # The final equation prints each number explicitly as requested.
    print(f"   - Total minimum stickers = {corner_stickers_needed} (corners) + {edge_stickers_needed} (edges) = {total_stickers}")

solve_pll_stickers()