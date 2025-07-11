# A python script to count the colors of phosphorus allotropes.

# The main, well-characterized allotropes of phosphorus and their colors.
# Note: White phosphorus is often called yellow phosphorus.
phosphorus_colors = {
    "White (or Yellow) Phosphorus": "white/yellow",
    "Red Phosphorus": "red",
    "Violet Phosphorus": "violet",
    "Black Phosphorus": "black"
}

# Get the list of unique colors.
colors = list(phosphorus_colors.values())
count = len(colors)

print("The main pure allotropes of phosphorus exhibit several distinct colors:")
for allotrope, color in phosphorus_colors.items():
    print(f"- {allotrope} is {color}.")

print("\nTo find the total, we can count each distinct color:")

# Create a simple sum equation (e.g., 1 + 1 + 1 + 1) to show the count.
# This represents counting each of the colors listed.
equation_numbers = ['1'] * count
equation_string = " + ".join(equation_numbers)

print(f"Equation: {equation_string} = {count}")

print(f"\nTherefore, {count} distinct colors can be observed in pure allotropes of phosphorus.")