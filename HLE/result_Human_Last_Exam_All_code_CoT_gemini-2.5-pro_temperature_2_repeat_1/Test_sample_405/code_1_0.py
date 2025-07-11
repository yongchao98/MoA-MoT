def count_phosphorus_colors():
    """
    Calculates and prints the number of colors observed in pure allotropes of phosphorus.
    """
    # The common allotropes of phosphorus and their colors:
    # 1. White Phosphorus: Can appear as a waxy white or yellowish solid.
    # 2. Red Phosphorus: Has a red/dark red color.
    # 3. Violet Phosphorus (Hittorf's phosphorus): Appears violet or purple.
    # 4. Black Phosphorus: Has a black, flaky appearance, similar to graphite.
    
    colors = ["White", "Yellow", "Red", "Violet", "Black"]
    number_of_colors = len(colors)

    print("The distinct colors of pure phosphorus allotropes are:")
    
    # Build the equation string as requested
    equation_parts = []
    for color in colors:
        print(f"- {color}")
        equation_parts.append(f"1")

    equation_str = " + ".join(equation_parts)

    print(f"\nSumming the number of distinct colors:")
    print(f"{equation_str} = {number_of_colors}")
    print(f"\nThere are {number_of_colors} different colors that can be observed in pure allotropes of phosphorus.")

count_phosphorus_colors()