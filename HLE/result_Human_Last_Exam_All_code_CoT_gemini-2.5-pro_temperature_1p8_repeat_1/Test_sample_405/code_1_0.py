# The main allotropes of phosphorus each exhibit a distinct color.
# This script will list them and provide a count.

# A dictionary mapping each allotrope to its characteristic color.
allotropes_and_colors = {
    "White Phosphorus": "White (or yellowish)",
    "Red Phosphorus": "Red",
    "Violet Phosphorus": "Violet",
    "Black Phosphorus": "Black",
    "Blue Phosphorus": "Blue" # A 2D allotrope, similar to graphene
}

print("The main allotropes of phosphorus and their corresponding colors are:")
for allotrope, color in allotropes_and_colors.items():
    print(f"- {color}")

# The total number of colors is the number of allotropes listed.
count = len(allotropes_and_colors)

# Building the equation string as 1 + 1 + ... for each color
equation_numbers = ["1"] * count
equation_string = " + ".join(equation_numbers)

print("\nTo find the total number of colors, we sum one for each distinct allotrope:")
print(f"{equation_string} = {count}")