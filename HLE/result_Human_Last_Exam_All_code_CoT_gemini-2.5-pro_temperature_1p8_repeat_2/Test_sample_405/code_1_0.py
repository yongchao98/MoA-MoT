# A list of the allotropes of phosphorus with distinct colors.
# While white phosphorus can appear yellowish, it's counted as one color.
# Violet (Hittorf's) phosphorus is distinct from the more common red phosphorus.
# Black phosphorus has two forms (orthorhombic and rhombohedral) but both are black.
phosphorus_allotropes = ["White", "Red", "Violet", "Black"]

# The number of colors is the number of items in our list.
total_colors = len(phosphorus_allotropes)

print("The number of distinct colors for the main allotropes of phosphorus can be found by summing them up.")
print("The allotropes are: White, Red, Violet, and Black.")
print("The calculation is:")

# Build and print the equation string, e.g., "1 + 1 + 1 + 1 = 4"
equation_numbers = ["1" for _ in phosphorus_allotropes]
equation_str = " + ".join(equation_numbers)
print(f"{equation_str} = {total_colors}")