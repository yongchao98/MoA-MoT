def count_phosphorus_colors():
    """
    This function identifies and counts the distinct colors of pure phosphorus allotropes.
    """
    # The primary colors of well-established phosphorus allotropes
    # - White Phosphorus is a white/colorless waxy solid when pure.
    # - Red Phosphorus is an amorphous polymer with a red color.
    # - Violet (or Hittorf's) Phosphorus has a deep violet crystalline structure.
    # - Black Phosphorus is the most stable form and is black with a metallic sheen.
    # - Blue Phosphorus is a more recently synthesized 2D allotrope.
    colors = [
        "White",
        "Red",
        "Violet",
        "Black",
        "Blue"
    ]

    print("The distinct colors observed in pure allotropes of phosphorus are:")
    for color in colors:
        print(f"- {color}")

    total_colors = len(colors)
    
    # Constructing the equation string as requested
    equation_parts = ["1" for color in colors]
    equation_string = " + ".join(equation_parts)
    
    print("\nCalculating the total number of colors:")
    print(f"1 (White) + 1 (Red) + 1 (Violet) + 1 (Black) + 1 (Blue) = {total_colors}")

count_phosphorus_colors()
<<<5>>>