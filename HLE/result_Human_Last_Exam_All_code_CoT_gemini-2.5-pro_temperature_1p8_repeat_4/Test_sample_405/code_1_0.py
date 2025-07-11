# The major allotropes of phosphorus and their distinct colors.
colors_of_phosphorus = {
    "White Phosphorus": "White/Yellow",
    "Red Phosphorus": "Red",
    "Violet Phosphorus": "Violet",
    "Black Phosphorus": "Black",
    "Blue Phosphorus": "Blue"
}

# The number of distinct colors is the number of allotropes listed.
count = len(colors_of_phosphorus)

# Build the equation string
equation_parts = []
for allotrope in colors_of_phosphorus:
    equation_parts.append(f"1 ({allotrope})")

equation_string = " + ".join(equation_parts)

# Print the final result
print("The calculation for the number of observable colors in phosphorus allotropes is:")
print(f"{equation_string} = {count}")
<<<5>>>