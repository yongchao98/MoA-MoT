def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    observed in the pure allotropes of phosphorus.
    """
    # A dictionary mapping phosphorus allotropes or categories to their observed colors.
    # Some allotropes can exhibit more than one color.
    allotrope_colors = {
        "White Phosphorus": ["White", "Yellow"],
        "Red Phosphorus": ["Red", "Scarlet"],
        "Violet Phosphorus": ["Violet"],
        "Black Phosphorus": ["Black"],
        "2D Phosphorus": ["Blue", "Green"] # Newer, synthesized 2D forms
    }

    print("The allotropes of phosphorus and their associated colors are:")
    
    # Collect all unique colors from the dictionary values
    unique_colors = set()
    for allotrope, colors in allotrope_colors.items():
        print(f"- {allotrope}: {', '.join(colors)}")
        for color in colors:
            unique_colors.add(color)

    # Sort for consistent ordering
    sorted_colors = sorted(list(unique_colors))
    
    total_colors = len(sorted_colors)
    
    print(f"\nThe unique colors are: {', '.join(sorted_colors)}.")
    print(f"This gives a total of {total_colors} observable colors.")
    
    # Building the equation string as requested
    # It represents counting each unique color once.
    equation_parts = ["1" for _ in sorted_colors]
    equation_str = " + ".join(equation_parts)
    
    print("\nThe final calculation is:")
    print(f"{equation_str} = {total_colors}")

# Execute the function
count_phosphorus_colors()