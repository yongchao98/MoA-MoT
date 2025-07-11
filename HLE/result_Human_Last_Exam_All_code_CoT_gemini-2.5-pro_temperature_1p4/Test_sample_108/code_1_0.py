import math

# Define the lengths of the robot arm segments in cm
l1 = 40.0  # Shoulder to Elbow
l2 = 28.0  # Elbow to Wrist
l3 = 15.0  # Wrist to Hand
l4 = 10.0  # Hand (finger)

# Define the clearance constraint at each joint in cm
# This is interpreted as the vertical separation for each folding segment
y_sep = 3.5

# Calculate the horizontal component of each folded segment's vector
# x = sqrt(L^2 - y^2)
x2 = math.sqrt(l2**2 - y_sep**2)
x3 = math.sqrt(l3**2 - y_sep**2)
x4 = math.sqrt(l4**2 - y_sep**2)

# Calculate the final coordinates of the fingertip
# The arm starts by extending 40cm in the x-direction (l1),
# then folds back with segments l2, l3, and l4.
final_x = l1 - x2 - x3 - x4

# The final y-coordinate is the sum of the vertical separations
final_y = y_sep + y_sep + y_sep

# Calculate the final distance from the shoulder (origin) to the fingertip
distance = math.sqrt(final_x**2 + final_y**2)

# Print the final equation with all numbers
print(f"To find the final position, we model the arm's folding in a plane.")
print(f"The first segment L1 extends along the x-axis by {l1} cm.")
print(f"The next three segments fold back, each with a vertical separation of {y_sep} cm to avoid collision.")
print(f"The horizontal projections (x_i) of these folded segments are calculated as x_i = sqrt(L_i^2 - y_sep^2).")
print(f"x2 = sqrt({l2}^2 - {y_sep}^2) = {x2:.4f} cm")
print(f"x3 = sqrt({l3}^2 - {y_sep}^2) = {x3:.4f} cm")
print(f"x4 = sqrt({l4}^2 - {y_sep}^2) = {x4:.4f} cm")
print(f"The final X coordinate is: Fx = {l1} - {x2:.4f} - {x3:.4f} - {x4:.4f} = {final_x:.4f} cm")
print(f"The final Y coordinate is: Fy = {y_sep} + {y_sep} + {y_sep} = {final_y:.4f} cm")
print("\nThe final distance from the shoulder (origin) to the finger is the magnitude of this position vector:")
print(f"Distance = sqrt(Fx^2 + Fy^2)")
print(f"Distance = sqrt({final_x:.4f}^2 + {final_y:.4f}^2)")
print(f"Distance = {distance:.4f} cm")

# The result is very close to one of the options.
print("\nThis calculated distance is approximately 15.75 cm.")
