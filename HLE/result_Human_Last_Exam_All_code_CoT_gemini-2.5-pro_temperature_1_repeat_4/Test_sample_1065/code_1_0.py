import math

# The final formulas for the horizontal (x) and vertical (y) coordinates of the center of mass
# are expressed as multiples of the radius R.
# x_cm / R = 4 / (pi + 2)
# y_cm / R = 1 / (pi + 2)

# The numbers that form these equations are 4, 1, and 2.
num_x_numerator = 4
num_y_numerator = 1
num_denominator_add = 2

# Calculate the numerical values of the coordinates (as raw numbers, i.e., coefficients of R)
pi = math.pi
horizontal_coordinate = num_x_numerator / (pi + num_denominator_add)
vertical_coordinate = num_y_numerator / (pi + num_denominator_add)

# Print the raw numbers for the horizontal and vertical coordinates, separated by a comma.
print(f"{horizontal_coordinate},{vertical_coordinate}")
