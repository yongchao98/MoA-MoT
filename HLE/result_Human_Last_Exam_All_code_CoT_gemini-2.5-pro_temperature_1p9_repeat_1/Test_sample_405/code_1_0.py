def count_phosphorus_colors():
    """
    This function identifies and counts the number of observable colors 
    in the pure allotropes of phosphorus.
    """
    
    # The distinct observable colors are derived from the main allotropes:
    # - White Phosphorus (can appear white or yellow)
    # - Red Phosphorus (is red)
    # - Violet Phosphorus (is violet)
    # - Black Phosphorus (is black)
    colors = ["White", "Yellow", "Red", "Violet", "Black"]
    
    # Get the total count of the colors.
    total_count = len(colors)
    
    print(f"The distinct colors that can be observed in pure allotropes of phosphorus are: {', '.join(colors)}.")
    
    # Create the left side of the equation, representing each color as '1'.
    equation_left_side = " + ".join(['1'] * total_count)
    
    # Print the final equation showing each number being added to get the total.
    print("The summation to find the total number of colors is:")
    print(f"{equation_left_side} = {total_count}")

count_phosphorus_colors()