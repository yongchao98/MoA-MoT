def count_phosphorus_colors():
    """
    Identifies and counts the distinct colors of the main allotropes of phosphorus.
    """
    # The main allotropes of phosphorus and their corresponding colors are:
    # 1. White Phosphorus (can also appear yellow)
    # 2. Red Phosphorus
    # 3. Violet Phosphorus (Hittorf's phosphorus)
    # 4. Black Phosphorus
    colors = ["White", "Red", "Violet", "Black"]

    # Calculate the total number of colors.
    total_colors = len(colors)

    print("The colors observed in the main pure allotropes of phosphorus are:")
    for i, color in enumerate(colors):
        print(f"{i+1}. {color}")

    # To satisfy the "output each number in the final equation" requirement,
    # we can represent the count as a sum of 1 for each color found.
    equation_numbers = ["1"] * len(colors)
    equation_str = " + ".join(equation_numbers)

    print(f"\nThe final count can be represented as the sum of each observed color:")
    print(f"{equation_str} = {total_colors}")

    print(f"\nIn total, {total_colors} distinct colors can be observed in pure allotropes of phosphorus.")

count_phosphorus_colors()