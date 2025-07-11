def count_phosphorus_colors():
    """
    Identifies and counts the distinct colors of phosphorus allotropes.
    """
    # A dictionary mapping phosphorus forms to their observed colors.
    # White phosphorus is included twice for its pure color and the color it turns upon exposure to light.
    allotropes_to_colors = {
        "White Phosphorus": "white",
        "White Phosphorus (aged/exposed)": "yellow",
        "Red Phosphorus": "red",
        "Violet Phosphorus (Hittorf's)": "violet",
        "Black Phosphorus": "black"
    }

    # Extract the unique colors from the dictionary values
    distinct_colors = sorted(list(set(allotropes_to_colors.values())))
    
    # Count the number of unique colors
    num_colors = len(distinct_colors)

    # Build the equation string to show the calculation
    # We represent each color as '1' in a sum.
    equation_parts = []
    for color in distinct_colors:
        equation_parts.append(f"1 ({color})")
    
    equation_str = " + ".join(equation_parts)
    
    print("The distinct colors of phosphorus allotropes are:")
    for color in distinct_colors:
        print(f"- {color.capitalize()}")
    
    print("\nThe equation to find the total number of colors is:")
    print(f"{equation_str} = {num_colors}")

count_phosphorus_colors()