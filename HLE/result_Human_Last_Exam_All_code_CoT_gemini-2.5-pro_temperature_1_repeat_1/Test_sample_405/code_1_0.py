def solve_phosphorus_colors():
    """
    This function identifies and counts the number of observable colors
    in the pure allotropes of phosphorus.
    """
    # The main allotropes of phosphorus have distinct colors.
    # We list them here.
    colors = {
        "White/Yellow": "White Phosphorus",
        "Red": "Red Phosphorus",
        "Violet": "Violet Phosphorus",
        "Black": "Black Phosphorus"
    }

    # The number of colors is the number of items in our dictionary.
    num_colors = len(colors)

    print("The distinct colors observed in pure phosphorus allotropes are:")
    for color, allotrope in colors.items():
        print(f"- {color} (from {allotrope})")

    print("\nTo find the total number of colors, we can represent each color as '1' and sum them up.")

    # Building and printing the equation as requested.
    equation_parts = ["1" for color in colors]
    equation_str = " + ".join(equation_parts)
    print(f"The equation is: {equation_str} = {num_colors}")

    print(f"\nThus, there are {num_colors} distinct colors that can be observed.")

solve_phosphorus_colors()