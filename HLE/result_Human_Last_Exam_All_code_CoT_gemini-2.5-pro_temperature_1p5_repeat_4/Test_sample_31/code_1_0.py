def solve_pll_sticker_count():
    """
    Explains the reasoning for the number of stickers needed to identify a PLL case.
    """
    
    # The 4 corner pieces of the top layer.
    corners = ["Front-Up-Left", "Front-Up-Right", "Right-Up-Back", "Back-Up-Left"]
    
    # Each corner piece has 2 non-top-facing stickers.
    stickers_per_corner = 2
    
    # The total number of non-top-facing corner stickers.
    # This is the information required to know the permutation of all corners.
    total_corner_stickers = len(corners) * stickers_per_corner
    
    print("To fully identify a PLL case, we need to distinguish it from all 20 other possibilities.")
    print("A robust method for this is to first determine the permutation of the 4 top-layer corner pieces.")
    print("\nEach of the 4 corner pieces has 2 side-facing stickers.")
    print(f"There are {len(corners)} corners and {stickers_per_corner} side stickers per corner.")
    print("To know where every corner is and how it's oriented, we must see all of these stickers.")
    
    print(f"\nThe calculation is: {len(corners)} corners * {stickers_per_corner} stickers/corner = {total_corner_stickers} stickers.")
    
    print("\nSeeing these 8 stickers allows you to identify patterns like 'headlights' on all four faces, which is the worst-case requirement for identifying the corner permutation.")
    print("Once the corner permutation is known, the PLL case is either uniquely identified or narrowed down to a very small group that can be distinguished by seeing just one or two more edge stickers.")
    print("\nTherefore, in the most challenging cases, 8 non-top-facing stickers must be seen to begin the identification.")
    
    # Printing the final answer in the requested format
    # The numbers in the equation are explicitly referenced.
    num_corners = len(corners)
    print(f"\nFinal Answer Equation: {num_corners} * {stickers_per_corner} = {total_corner_stickers}")


solve_pll_sticker_count()

# The final answer as a raw value
# <<<8>>>