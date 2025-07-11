def count_phosphorus_colors():
    """
    Identifies and counts the distinct colors of phosphorus allotropes.
    """
    # A dictionary mapping each distinct color to the allotrope(s) that exhibit it.
    allotropes = {
        "White": "White Phosphorus",
        "Yellow": "White Phosphorus (when exposed to light)",
        "Red": "Red Phosphorus",
        "Violet": "Violet (Hittorf's) Phosphorus",
        "Black": "Black Phosphorus"
    }

    print("The main colors observed in the pure allotropes of phosphorus are:")
    # Print each color and its corresponding allotrope.
    for color, name in allotropes.items():
        print(f"- {color} (from {name})")

    # Get the total number of distinct colors.
    num_colors = len(allotropes)
    
    # Create the equation string as a sum of 1s.
    equation_parts = ["1"] * num_colors
    equation_str = " + ".join(equation_parts)

    print("\nCalculating the total number of distinct colors:")
    print(f"{equation_str} = {num_colors}")

# Run the function to display the result.
count_phosphorus_colors()