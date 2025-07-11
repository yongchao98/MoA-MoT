def count_phosphorus_colors():
    """
    This function calculates and prints the number of observable colors
    in the pure allotropes of phosphorus.
    """
    # A dictionary mapping phosphorus allotropes to their observed colors.
    # White phosphorus can be white or yellow.
    allotropes_and_colors = {
        "White Phosphorus": ["white", "yellow"],
        "Red Phosphorus": ["red"],
        "Violet Phosphorus": ["violet"],
        "Black Phosphorus": ["black"]
    }

    print("The main allotropes of phosphorus and their colors are:")
    
    all_colors = []
    counts_per_allotrope = []

    # Iterate through the dictionary to display the information and collect counts.
    for allotrope, colors in allotropes_and_colors.items():
        print(f"- {allotrope}: {len(colors)} color(s) ({', '.join(colors)})")
        counts_per_allotrope.append(str(len(colors)))
        all_colors.extend(colors)

    # Calculate the total number of colors.
    total_colors = len(all_colors)

    # Create the equation string as requested.
    equation_string = " + ".join(counts_per_allotrope)

    print("\nThe total number of colors is the sum of the number of colors for each allotrope:")
    print(f"{equation_string} = {total_colors}")

# Execute the function.
count_phosphorus_colors()
<<<5>>>