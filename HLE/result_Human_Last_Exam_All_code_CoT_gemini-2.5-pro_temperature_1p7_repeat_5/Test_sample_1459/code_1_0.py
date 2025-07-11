import math

# This script calculates the Gromov-Hausdorff distance between the interval [0,1]
# and the unit circle with their intrinsic metrics.

# Define the parameters of the spaces
# L is the length of the interval.
L = 1.0
# R is the radius of the unit circle.
R = 1.0
# C is the circumference of the unit circle.
C = 2 * math.pi * R

# According to a known formula, the Gromov-Hausdorff distance d_GH is:
# d_GH = (L + C/2 - d) / 2
# where d is the length of a chord connecting two points at arc distance L.

# On a circle of radius R, an arc of length L subtends an angle theta = L/R.
theta = L / R

# The length of the chord d is given by d = 2 * R * sin(theta / 2).
d = 2 * R * math.sin(theta / 2)

# Substitute the values into the formula for d_GH.
d_gh = (L + C / 2 - d) / 2

# We need to print the full equation with the numerical values.
L_val_str = str(L)
pi_val_str = str(math.pi)
d_val_str = str(d)
numerator_val_str = str(L + math.pi - d)
final_result_str = str(d_gh)

print("The final equation for the Gromov-Hausdorff distance is:")
# The symbolic formula is D = (L + pi - 2*sin(L/(2*R))) / 2
# Substituting L=1, R=1, we get D = (1 + pi - 2*sin(1/2)) / 2
# The code below prints the equation with evaluated numbers.
print(f"({L_val_str} + {pi_val_str} - {d_val_str}) / 2 = {final_result_str}")
