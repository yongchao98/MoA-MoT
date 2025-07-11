def count_phosphorus_colors():
    """
    This function calculates the number of distinct colors observed
    in the pure allotropes of phosphorus.
    """
    # Dictionary mapping phosphorus allotropes to their observed colors.
    # White phosphorus is often called yellow phosphorus as it turns yellow in light.
    allotropes = {
        "White/Yellow Phosphorus": ["white", "yellow"],
        "Red Phosphorus": ["red"],
        "Violet Phosphorus": ["violet"],
        "Black Phosphorus": ["black"]
    }

    print("Calculating the number of observable colors in pure phosphorus allotropes.")
    print("-" * 60)

    all_colors = []
    equation_parts = []

    # Iterate through the allotropes to gather colors and explain the count.
    for allotrope, colors in allotropes.items():
        count = len(colors)
        print(f"{allotrope}: Contributes {count} color(s) ({', '.join(colors)})")
        all_colors.extend(colors)
        equation_parts.append(str(count))

    # Using a set to find the number of unique colors.
    unique_colors = set(all_colors)
    total_unique_colors = len(unique_colors)

    # Building and printing the final equation as requested.
    final_equation = " + ".join(equation_parts)
    print("-" * 60)
    print("The final calculation is based on the number of colors from each category:")
    print(f"{final_equation} = {total_unique_colors}")
    print("-" * 60)
    print(f"There are a total of {total_unique_colors} distinct colors observed in the pure allotropes of phosphorus.")

if __name__ == "__main__":
    count_phosphorus_colors()