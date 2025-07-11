def count_phosphorus_colors():
    """
    This function identifies and counts the distinct colors of phosphorus allotropes.
    """
    # Allotropes of phosphorus exhibit several distinct colors.
    # 1. White Phosphorus: Can be white or yellow.
    # 2. Red Phosphorus: Is red.
    # 3. Violet Phosphorus: Is violet.
    # 4. Black Phosphorus: Is black.
    colors = ["White", "Yellow", "Red", "Violet", "Black"]
    
    # Calculate the total number of distinct colors.
    total_colors = len(colors)
    
    # Create the equation string.
    # Each color is counted as one item in our list.
    equation_parts = []
    for color in colors:
        equation_parts.append(f"1 ({color})")
    
    equation_string = " + ".join(equation_parts)
    
    print(f"The number of observable colors in pure phosphorus allotropes can be calculated as follows:")
    print(f"{equation_string} = {total_colors}")
    print(f"\nTherefore, there are {total_colors} distinct colors that can be observed.")

count_phosphorus_colors()