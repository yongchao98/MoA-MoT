import math

# We will calculate a proven lower bound for the radius of a circle
# containing 14 unit circles.
# The formula is R_min = 1 + 1 / sin(pi/9).

# First, define the angle in radians.
angle_in_radians = math.pi / 9

# Calculate the sine of the angle.
sin_of_angle = math.sin(angle_in_radians)

# Now, apply the formula to find the radius.
radius = 1 + (1 / sin_of_angle)

# As requested, here is the final equation with each number shown.
print("The formula for the lower bound of the radius (R) is: R = 1 + 1 / sin(pi/9)")
print("Plugging in the values, the equation is:")
print(f"R = 1 + 1 / {sin_of_angle} = {radius}")

# Finally, we provide the answer rounded to 4 significant digits.
# The .4g format specifier rounds to 4 significant figures.
radius_rounded = float(f"{radius:.4g}")

print(f"\nThe resulting radius, up to 4 significant digits, is {radius_rounded}.")