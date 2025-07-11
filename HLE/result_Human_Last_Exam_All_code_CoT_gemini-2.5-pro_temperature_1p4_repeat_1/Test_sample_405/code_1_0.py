def count_phosphorus_colors():
    """
    This function identifies the colors of the main allotropes of phosphorus
    and calculates the total number of distinct colors.
    """
    # Step 1 & 2: Define the main allotropes and their distinct colors.
    # While white phosphorus turns yellow on exposure to light, it represents one color type.
    allotropes = {
        "White (or Yellow) Phosphorus": 1,
        "Red Phosphorus": 1,
        "Violet Phosphorus": 1,
        "Black Phosphorus": 1
    }

    print("The main allotropes of phosphorus exhibit the following colors:")
    for color in allotropes.keys():
        print(f"- {color}")

    # Step 3: Count the number of unique colors.
    num_colors = len(allotropes)
    
    # Step 4: Create and print the equation as requested.
    # We create a list of '1's, one for each color.
    summands = ['1'] * num_colors
    
    # We join the '1's with a '+' sign to form the left side of the equation.
    equation_str = " + ".join(summands)

    print("\nCounting each distinct color, the total is:")
    print(f"{equation_str} = {num_colors}")

count_phosphorus_colors()