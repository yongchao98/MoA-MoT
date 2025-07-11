def count_phosphorus_colors():
    """
    This function identifies the colors of the main allotropes of phosphorus
    and counts the number of distinct colors.
    """
    # Step 1 & 2: Define the main allotropes and their corresponding colors.
    # We focus on the four major, well-established allotropes.
    # - White Phosphorus is a molecular solid (P4), appearing white or colorless when pure.
    # - Red Phosphorus is a polymeric form, appearing red.
    # - Violet Phosphorus (Hittorf's Phosphorus) is another polymeric form, appearing violet.
    # - Black Phosphorus is the most stable form, layered like graphite, appearing black.
    allotropes_to_colors = {
        "White Phosphorus": "white",
        "Red Phosphorus": "red",
        "Violet Phosphorus": "violet",
        "Black Phosphorus": "black"
    }

    # Step 3: Get the list of unique colors.
    unique_colors = sorted(list(set(allotropes_to_colors.values())))
    
    # Step 4: Display the findings and the final count.
    # We will list each color found and then sum them for the total.
    print("The distinct colors observed in the main pure allotropes of phosphorus are:")
    for color in unique_colors:
        print(f"- {color.capitalize()}")

    print("\nTo find the total number of colors, we sum the count for each unique color:")
    
    # Building the "equation" as requested
    equation_parts = [f"1 ({color})" for color in unique_colors]
    equation_string = " + ".join(equation_parts)
    
    total = len(unique_colors)
    
    print(f"{equation_string} = {total}")

if __name__ == "__main__":
    count_phosphorus_colors()