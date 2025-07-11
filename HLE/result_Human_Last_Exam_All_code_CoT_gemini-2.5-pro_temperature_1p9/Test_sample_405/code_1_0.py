def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    observed in the pure allotropes of phosphorus.
    """
    # Dictionary mapping the allotrope name to its characteristic color
    allotropes = {
        "White Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black"
    }

    print("The main allotropes of phosphorus and their colors are:")
    for allotrope, color in allotropes.items():
        print(f"- {allotrope}: {color}")

    count = len(allotropes)
    
    # Building the equation string as requested
    equation_parts = ["1" for _ in allotropes]
    equation_str = " + ".join(equation_parts)
    
    print("\nTo find the total number of colors, we can set up the following equation:")
    print(f"{equation_str} = {count}")

count_phosphorus_colors()