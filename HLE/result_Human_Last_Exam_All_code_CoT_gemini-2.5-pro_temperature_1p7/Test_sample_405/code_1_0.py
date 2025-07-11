def count_phosphorus_colors():
    """
    Identifies and counts the unique colors of phosphorus allotropes.
    """
    # A dictionary mapping phosphorus allotropes to their observed colors.
    # White phosphorus is often considered white/yellow.
    # Violet phosphorus is also known as Hittorf's phosphorus.
    # Blue phosphorus is a more recently discovered 2D allotrope.
    allotropes = {
        "White (or Yellow) Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black",
        "Blue Phosphorus": "blue"
    }

    print("The primary colors observed in pure allotropes of phosphorus are:")
    
    # Extract the unique colors
    colors = list(allotropes.values())
    
    # Print each color found
    for color in colors:
        print(f"- {color.capitalize()}")

    # Build the equation string as requested
    num_colors = len(colors)
    equation_parts = ["1"] * num_colors
    equation_str = " + ".join(equation_parts)

    print(f"\nCounting each unique color gives the equation:")
    print(f"{equation_str} = {num_colors}")

if __name__ == "__main__":
    count_phosphorus_colors()