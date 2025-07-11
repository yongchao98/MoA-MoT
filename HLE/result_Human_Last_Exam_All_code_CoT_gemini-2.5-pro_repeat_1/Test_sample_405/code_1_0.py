def count_phosphorus_colors():
    """
    This script identifies the colors of the common allotropes of phosphorus
    and calculates the total number of distinct colors.
    """
    
    # A dictionary mapping each common allotrope to its distinct color.
    allotropes = {
        "White Phosphorus": "White (or yellowish)",
        "Red Phosphorus": "Red",
        "Violet Phosphorus": "Violet",
        "Black Phosphorus": "Black"
    }
    
    print("The colors of the pure allotropes of phosphorus are:")
    for allotrope, color in allotropes.items():
        print(f"- {color} (from {allotrope})")
        
    num_colors = len(allotropes)
    
    print(f"\nThere are {num_colors} distinct colors that can be observed.")
    
    # Building and printing the final equation as requested.
    equation_parts = ["1"] * num_colors
    equation_string = " + ".join(equation_parts)
    
    print("The final equation for the total number of colors is:")
    print(f"{equation_string} = {num_colors}")

count_phosphorus_colors()