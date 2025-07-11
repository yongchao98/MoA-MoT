def count_phosphorus_colors():
    """
    This function identifies the colors of phosphorus allotropes,
    counts the unique ones, and prints the result in an equation format.
    """
    # A dictionary mapping phosphorus allotropes to their observed colors.
    # Note: Scarlet phosphorus is a form of red, and white can appear yellow,
    # but these are the fundamental distinct colors of the pure allotropes.
    allotrope_colors = {
        "White Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black",
        "Blue Phosphorus": "blue"
    }

    # Extract the color values and find the unique set of colors
    unique_colors = sorted(list(set(allotrope_colors.values())))

    # Count the number of unique colors
    num_colors = len(unique_colors)

    print("The distinct colors observed in pure allotropes of phosphorus are:")
    for color in unique_colors:
        print(f"- {color.capitalize()}")

    print("\nTo find the total, we sum the count for each unique color.")

    # Build and print the equation as requested
    # Each unique color contributes 1 to the total count.
    equation_numbers = ["1"] * num_colors
    equation_string = " + ".join(equation_numbers)

    print(f"Final Equation: {equation_string} = {num_colors}")

# Execute the function
count_phosphorus_colors()