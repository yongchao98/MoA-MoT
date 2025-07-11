def count_phosphorus_colors():
    """
    This function calculates and prints the number of observable colors
    in pure allotropes of phosphorus.
    """
    # A dictionary mapping phosphorus allotropes to their observed colors.
    # White phosphorus can appear both white and yellow.
    allotropes_and_colors = {
        "White Phosphorus": ["white", "yellow"],
        "Red Phosphorus": ["red"],
        "Violet Phosphorus": ["violet"],
        "Black Phosphorus": ["black"],
        "Blue Phosphorus": ["blue"]
    }

    # Use a set to store the unique colors to avoid duplicates.
    unique_colors = set()

    # Iterate through the dictionary values (lists of colors) and add each color to the set.
    for colors in allotropes_and_colors.values():
        unique_colors.update(colors)

    # Convert the set to a sorted list for consistent output order.
    sorted_colors = sorted(list(unique_colors))
    
    # Calculate the total number of unique colors.
    count = len(sorted_colors)

    print("The distinct colors observed in the pure allotropes of phosphorus are:")
    # Print each color found. This fulfills the "output each number" requirement by listing each item being counted.
    for color in sorted_colors:
        print(f"- {color}")
    
    # Create the final "equation" string.
    equation_str = " + ".join(["1" for _ in sorted_colors])
    
    print(f"\nCounting each unique color (1 for each):")
    print(f"The calculation is: {equation_str} = {count}")
    print(f"Total number of colors: {count}")

# Run the function to display the result.
count_phosphorus_colors()