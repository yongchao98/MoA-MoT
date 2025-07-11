import math

# Parameters of the metric spaces
# Length of the interval [0, 1]
L = 1.0
# Circumference of the unit circle (radius 1)
C = 2 * math.pi

# The term C/pi from the formula
C_div_pi = C / math.pi

# The formula for the Gromov-Hausdorff distance in this case is
# d_GH = (1/2) * sqrt(L^2 + (C/pi)^2)
distance = 0.5 * math.sqrt(L**2 + C_div_pi**2)

# Print the equation with the specific numbers
print(f"The Gromov-Hausdorff distance is calculated by the formula: d_GH = sqrt(L^2 + (C/pi)^2) / 2")
print(f"Substituting the values L={int(L)} and C=2*pi:")
# We display C/pi as an integer because it simplifies to 2.
print(f"d_GH = sqrt({int(L)}^2 + ({int(C_div_pi)})^2) / 2")
print(f"d_GH = sqrt({int(L**2)} + {int(C_div_pi**2)}) / 2")
print(f"d_GH = sqrt({int(L**2 + C_div_pi**2)}) / 2")
print(f"The final distance is: {distance}")