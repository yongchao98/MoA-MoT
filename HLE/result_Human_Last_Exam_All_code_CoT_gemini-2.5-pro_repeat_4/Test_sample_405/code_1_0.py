def count_phosphorus_colors():
    """
    Identifies and counts the number of distinct colors of phosphorus allotropes.
    """
    # A dictionary mapping phosphorus allotropes/forms to their observed colors.
    # White phosphorus can appear as white or yellow.
    allotrope_colors = {
        "White Phosphorus": "White",
        "Yellow Phosphorus": "Yellow",
        "Red Phosphorus": "Red",
        "Violet Phosphorus": "Violet",
        "Black Phosphorus": "Black"
    }

    # Extract the unique color names from the dictionary values
    unique_colors = sorted(list(set(allotrope_colors.values())))
    
    # Get the total count of unique colors
    color_count = len(unique_colors)

    print("The distinct colors observed in pure allotropes of phosphorus are:")
    for color in unique_colors:
        print(f"- {color}")

    print("\nThe final count is derived from the sum of each unique color found.")
    
    # Create and print the equation as requested
    equation_parts = ["1"] * color_count
    equation_str = " + ".join(equation_parts)
    
    print(f"Equation: {equation_str} = {color_count}")

# Execute the function
count_phosphorus_colors()