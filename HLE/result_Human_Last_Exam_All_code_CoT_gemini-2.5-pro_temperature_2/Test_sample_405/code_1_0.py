def count_phosphorus_colors():
    """
    This function identifies the colors of phosphorus allotropes and calculates the total number.
    """
    # The main allotropes of phosphorus have distinct colors.
    # Let's list them to find the total count.
    colors_info = {
        "White phosphorus": "White (or yellowish)",
        "Red phosphorus": "Red",
        "Violet phosphorus": "Violet (Hittorf's phosphorus)",
        "Black phosphorus": "Black"
    }

    print("The distinct colors of pure phosphorus allotropes are:")
    for allotrope, color in colors_info.items():
        print(f"- {color}, from {allotrope}")

    num_colors = len(colors_info)

    # To fulfill the requirement of showing the numbers in an equation,
    # we will represent the count of each color as '1'.
    print("\nCalculating the total number of colors:")
    
    # Create the equation string like "1 + 1 + 1 + 1 = 4"
    equation_parts = ["1"] * num_colors
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {num_colors}")
    print(f"\nIn total, there are {num_colors} primary colors that can be observed in pure allotropes of phosphorus.")

count_phosphorus_colors()