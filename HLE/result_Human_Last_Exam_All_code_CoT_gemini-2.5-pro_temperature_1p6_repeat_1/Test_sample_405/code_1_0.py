def count_phosphorus_colors():
    """
    This function identifies the colors of common phosphorus allotropes,
    counts them, and prints the result as an equation.
    """
    # List of distinct colors for the well-established allotropes of phosphorus
    allotrope_colors = ["White", "Red", "Violet", "Black"]

    # Calculate the total number of colors
    num_colors = len(allotrope_colors)

    # Build the equation string for display
    equation_parts = []
    for color in allotrope_colors:
        equation_parts.append(f"1 (for {color})")
    
    equation_string = " + ".join(equation_parts)
    
    # Print the final equation with the total count
    print(f"The number of observable colors in pure phosphorus allotropes can be calculated as:")
    print(f"{equation_string} = {num_colors}")

count_phosphorus_colors()