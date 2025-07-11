def count_phosphorus_colors():
    """
    This function counts and displays the number of observable colors
    in the pure allotropes of phosphorus.
    """
    # Allotropes and their colors:
    # 1. White Phosphorus (can also appear Yellow)
    # 2. Red Phosphorus
    # 3. Violet Phosphorus (Hittorf's Phosphorus)
    # 4. Black Phosphorus
    # 5. Scarlet Phosphorus
    # 6. Blue Phosphorus (a 2D form)
    colors = ["White", "Yellow", "Red", "Violet", "Black", "Scarlet", "Blue"]
    
    num_colors = len(colors)
    
    print("The allotropes of phosphorus can be observed in several distinct colors.")
    print("The colors are: " + ", ".join(colors) + ".")
    
    # Building and printing the equation as requested
    equation_str = " + ".join(["1"] * num_colors)
    print(f"The total number of colors is derived from the following calculation:")
    print(f"{equation_str} = {num_colors}")

count_phosphorus_colors()