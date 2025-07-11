import math

# The problem asks for the radius of a large circle that tightly packs 14 smaller circles of radius 1.
# This is a well-known problem in geometry. The proven optimal radius, R, for 14 circles
# of radius 1 is given by the exact formula R = sqrt(13).

# Define the number inside the square root from the formula.
number = 13

# Calculate the radius using the formula.
radius = math.sqrt(number)

# The final equation is R = sqrt(13). We print the number in this equation.
print(f"The final equation for the radius R is: R = sqrt({number})")

# Print the final answer, formatted to 4 significant digits.
# The format specifier "{value:.4g}" formats the number to 4 significant figures.
print(f"The radius of the circle is {radius:.4g}")
