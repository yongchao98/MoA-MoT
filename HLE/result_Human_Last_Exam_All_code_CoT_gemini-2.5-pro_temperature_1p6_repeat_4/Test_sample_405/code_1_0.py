def count_phosphorus_colors():
    """
    This function identifies the colors of pure phosphorus allotropes,
    counts the number of distinct colors, and prints the result as an equation.
    """
    # A dictionary mapping major phosphorus allotropes to their colors.
    allotropes = {
        "White Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black"
    }

    print("The primary allotropes of phosphorus and their colors are:")
    for allotrope, color in allotropes.items():
        print(f"- {allotrope}: {color}")

    # Extract the unique colors from the dictionary values.
    # Using a set automatically handles uniqueness.
    unique_colors = set(allotropes.values())
    
    # Get the count of unique colors.
    num_colors = len(unique_colors)

    # Create the parts of the equation. Each unique color contributes '1' to the sum.
    equation_parts = ['1'] * num_colors
    
    # Join the parts with ' + ' to form the equation string.
    equation_str = " + ".join(equation_parts)

    print("\nTo find the total number of colors, we sum the count for each distinct color.")
    # The final print statement shows the equation with each number, as requested.
    print(f"The calculation is: {equation_str} = {num_colors}")
    print(f"\nTherefore, there are {num_colors} distinct colors that can be observed in pure allotropes of phosphorous.")

# Execute the function
count_phosphorus_colors()