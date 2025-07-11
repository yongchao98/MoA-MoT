def count_phosphorus_colors():
    """
    Identifies and counts the number of distinct colors of pure phosphorus allotropes.
    """
    # A dictionary mapping major allotropes of phosphorus to their observed colors.
    # While shades exist (e.g., white can turn yellow, red has a scarlet variant),
    # these are the four primary distinct color categories.
    allotropes = {
        "White Phosphorus": "White",
        "Red Phosphorus": "Red",
        "Violet Phosphorus": "Violet",
        "Black Phosphorus": "Black"
    }

    # Extract the unique colors from the dictionary values
    distinct_colors = sorted(list(set(allotropes.values())))

    # Get the total count of distinct colors
    color_count = len(distinct_colors)

    # Build the equation string as requested
    equation_parts = []
    for color in distinct_colors:
        equation_parts.append(f"1 ({color})")

    equation_string = " + ".join(equation_parts) + f" = {color_count}"

    print("The main colors observed in the pure allotropes of phosphorus are:")
    for color in distinct_colors:
        print(f"- {color}")

    print("\nTo find the total number of distinct colors, we can represent it as the following sum:")
    print(equation_string)

# Execute the function
count_phosphorus_colors()