def count_phosphorus_colors():
    """
    Identifies the colors of phosphorus allotropes and calculates the total number.
    """
    # Step 1 & 2: Define the main allotropes and their distinct colors.
    # White phosphorus can appear yellowish, but it's fundamentally one color category.
    allotropes = {
        "White Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black"
    }

    print("The primary allotropes of phosphorus exhibit the following colors:")
    for allotrope, color in allotropes.items():
        print(f"- {allotrope} is {color}.")

    # Step 3: Count the distinct colors.
    num_colors = len(allotropes)
    
    # Step 4: Display the final calculation as an equation.
    # The instruction requires showing each number in the equation.
    # We will represent each distinct color with the number '1'.
    equation_parts = ['1'] * num_colors
    equation_str = " + ".join(equation_parts)
    
    print("\nTo find the total number of colors, we sum one for each distinct color:")
    print(f"{equation_str} = {num_colors}")

count_phosphorus_colors()