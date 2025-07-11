def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    observed in the common allotropes of pure phosphorus.
    """
    
    # The distinct colors of the main phosphorus allotropes
    colors = {
        "White (or yellow)": 1,
        "Red": 1,
        "Violet": 1,
        "Black": 1
    }

    print("The distinct colors observed in the main allotropes of phosphorus are:")
    for color_name in colors.keys():
        print(f"- {color_name}")

    total_colors = len(colors)
    
    # Creating the equation string as requested
    equation_parts = ["1"] * total_colors
    equation_str = " + ".join(equation_parts)

    print(f"\nCounting these distinct colors gives the equation:")
    print(f"{equation_str} = {total_colors}")

count_phosphorus_colors()