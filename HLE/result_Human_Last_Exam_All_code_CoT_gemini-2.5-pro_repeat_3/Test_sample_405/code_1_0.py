def count_phosphorus_colors():
    """
    This function identifies the colors of phosphorus allotropes and counts them.
    """
    # A dictionary mapping major phosphorus allotropes to their colors.
    phosphorus_allotropes = {
        "White Phosphorus": "White",
        "Red Phosphorus": "Red",
        "Violet Phosphorus (Hittorf's)": "Violet",
        "Black Phosphorus": "Black"
    }

    print("The primary colors of pure phosphorus allotropes are:")
    # Using a set to get unique color values from the dictionary
    unique_colors = sorted(list(set(phosphorus_allotropes.values())))

    # Create a list of strings for the "equation" part
    equation_parts = []
    for color in unique_colors:
        print(f"- {color}")
        equation_parts.append(f"1 ({color})")

    # The final count is the number of unique colors
    total_count = len(unique_colors)

    # Join the parts to form a descriptive equation string
    final_equation_str = " + ".join(equation_parts)

    print("\nTo find the total number of distinct colors, we can form the equation:")
    print(f"{final_equation_str} = {total_count}")

# Execute the function
count_phosphorus_colors()