import math

# Segment lengths
L1 = 40  # Shoulder to Elbow
L2 = 28  # Elbow to Wrist
L3 = 15  # Wrist to Hand
L4 = 10  # Hand to Finger

# Collision constraint
min_separation = 1.0

# From the geometric analysis, we find the minimum sine of the joint angles
# to satisfy the collision constraint.
# sin(d_E) >= min_separation / L2
# sin(d_W) >= min_separation / L3
# d_H is unconstrained, so we set it to 0 to minimize distance.

sin_d_E = min_separation / L2
sin_d_W = min_separation / L3
sin_d_H = 0

# Calculate the cosine of the angles using sin^2 + cos^2 = 1
# cos = sqrt(1 - sin^2)
cos_d_E = math.sqrt(1 - sin_d_E**2)
cos_d_W = math.sqrt(1 - sin_d_W**2)
cos_d_H = math.sqrt(1 - sin_d_H**2) # This is 1

# Calculate the final coordinates of the finger tip
# The arm folds in a zig-zag pattern in the x-y plane
x_F = L1 - L2 * cos_d_E + L3 * cos_d_W - L4 * cos_d_H
y_F = L2 * sin_d_E + L3 * sin_d_W + L4 * sin_d_H

# Calculate the distance from the shoulder (origin) to the finger tip
distance = math.sqrt(x_F**2 + y_F**2)

# Output the equation and the final result
# The x-component of the final position is given by:
# x_F = L1 - L2*cos(d_E) + L3*cos(d_W) - L4*cos(d_H)
# We can express cos(d) as sqrt(1 - sin(d)^2)
# x_F = L1 - L2*sqrt(1 - (1/L2)^2) + L3*sqrt(1 - (1/L3)^2) - L4
# x_F = L1 - sqrt(L2^2 - 1) + sqrt(L3^2 - 1) - L4
x_eq_part1 = L1
x_eq_part2 = math.sqrt(L2**2 - 1)
x_eq_part3 = math.sqrt(L3**2 - 1)
x_eq_part4 = L4

# The y-component of the final position is given by:
# y_F = L2*sin(d_E) + L3*sin(d_W) + L4*sin(d_H)
# y_F = L2*(1/L2) + L3*(1/L3) + L4*0 = 1 + 1 + 0 = 2
y_eq_part1 = L2 * sin_d_E
y_eq_part2 = L3 * sin_d_W
y_eq_part3 = L4 * sin_d_H

print("The final distance is calculated using the Pythagorean theorem: distance = sqrt(x_F^2 + y_F^2)")
print(f"x_F = {x_eq_part1} - sqrt({L2}^2 - 1^2) + sqrt({L3}^2 - 1^2) - {x_eq_part4} = {x_F:.4f}")
print(f"y_F = {y_eq_part1:.0f} + {y_eq_part2:.0f} + {y_eq_part3:.0f} = {y_F:.4f}")
print(f"Distance = sqrt({x_F:.4f}^2 + {y_F:.4f}^2) = {distance:.2f} cm")