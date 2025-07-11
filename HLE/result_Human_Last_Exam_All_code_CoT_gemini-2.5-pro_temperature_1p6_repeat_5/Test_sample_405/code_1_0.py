# The main allotropes of phosphorus and their distinct colors
allotropes_with_colors = {
    "White Phosphorus": "White",
    "Red Phosphorus": "Red",
    "Black Phosphorus": "Black",
    "Violet Phosphorus": "Violet",
    "Blue Phosphorus": "Blue"
}

# Get the count of distinct colors
color_count = len(allotropes_with_colors)

print(f"There are {color_count} main observable colors in the allotropes of pure phosphorus.")
print("They are:")
for allotrope, color in allotropes_with_colors.items():
    print(f"- {color} (from {allotrope})")

# Build and print the final equation showing the count
equation_parts = [f"1 ({color})" for color in allotropes_with_colors.values()]
equation_string = " + ".join(equation_parts)
final_equation = f"\n{equation_string} = {color_count}"

print(final_equation)