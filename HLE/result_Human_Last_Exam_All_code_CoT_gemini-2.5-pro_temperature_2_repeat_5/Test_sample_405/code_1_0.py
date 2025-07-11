def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    observed in the common pure allotropes of phosphorus.
    """
    # The primary colors of the main allotropes of phosphorus.
    # 1. White Phosphorus (is often yellowish)
    # 2. Red Phosphorus
    # 3. Violet Phosphorus (Hittorf's Phosphorus)
    # 4. Black Phosphorus
    colors = [
        "White/Yellow",
        "Red",
        "Violet",
        "Black"
    ]
    
    # Count the number of colors
    num_colors = len(colors)
    
    # Create the equation string, showing '1' for each color found.
    equation_parts = ["1"] * num_colors
    equation_str = " + ".join(equation_parts)
    
    # Print the findings
    print(f"The distinct colors of the common pure allotropes of phosphorus are: {', '.join(colors)}.")
    print(f"Therefore, the total number of colors is {equation_str} = {num_colors}.")

count_phosphorus_colors()