# The main allotropes of phosphorus exhibit distinct colors.
# We will list these colors and then count them.

# List of colors for the main phosphorus allotropes.
colors = ["white", "red", "violet", "black"]

# Get the total number of colors.
num_colors = len(colors)

print(f"The colors of pure phosphorus allotropes are {', '.join(colors)}.")
print("This gives us a total count based on the following equation:")

# Create and print the equation as requested.
# "1" represents the count for each distinct color.
equation_parts = ["1"] * num_colors
equation_str = " + ".join(equation_parts)
print(f"{equation_str} = {num_colors}")

print(f"\nTherefore, there are {num_colors} distinct colors that can be observed in pure allotropes of phosphorus.")