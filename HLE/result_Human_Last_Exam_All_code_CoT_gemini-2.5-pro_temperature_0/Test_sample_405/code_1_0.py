def count_phosphorus_colors():
    """
    This function identifies and counts the colors of the main
    allotropes of phosphorus.
    """
    # A list of the distinct colors of phosphorus allotropes.
    colors = [
        "White (or yellow)",
        "Red",
        "Violet",
        "Black",
        "Blue"
    ]

    # Get the total number of colors.
    num_colors = len(colors)

    print("The distinct colors observed in the main allotropes of phosphorus are:")
    for color in colors:
        print(f"- {color}")

    # To represent the count as an equation, we show a sum of 1 for each color.
    equation_parts = ["1" for _ in colors]
    equation_str = " + ".join(equation_parts)

    print("\nTo find the total, we count each distinct color:")
    print(f"{equation_str} = {num_colors}")

    print(f"\nTherefore, {num_colors} distinct colors can be observed in pure allotropes of phosphorus.")

count_phosphorus_colors()