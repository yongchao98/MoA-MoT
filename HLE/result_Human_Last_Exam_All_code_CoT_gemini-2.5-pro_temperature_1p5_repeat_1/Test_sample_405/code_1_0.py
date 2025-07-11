def count_phosphorus_colors():
    """
    Identifies the colors of phosphorus allotropes and counts them.
    """
    # Step 1 & 2: Identify allotropes and their distinct colors.
    # Note: White phosphorus can turn yellow upon exposure to light,
    # but it's considered one color group (white/yellow).
    allotropes_colors = {
        "White/Yellow": 1,
        "Red": 1,
        "Violet": 1,
        "Black": 1
    }

    color_names = list(allotropes_colors.keys())
    total_colors = len(color_names)

    print("The main colors observed in pure allotropes of phosphorus are:")
    for color in color_names:
        print(f"- {color}")

    # Step 3 & 4: Create and print the equation representing the count.
    # Each '1' represents one of the distinct colors.
    equation_parts = [str(v) for v in allotropes_colors.values()]
    equation_string = " + ".join(equation_parts)

    print("\nTo find the total number of colors, we sum them up:")
    print(f"{equation_string} = {total_colors}")

count_phosphorus_colors()