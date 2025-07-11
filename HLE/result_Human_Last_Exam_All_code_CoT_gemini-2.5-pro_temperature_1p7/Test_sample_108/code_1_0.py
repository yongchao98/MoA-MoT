import math

# Define the segment lengths and clearance
L1 = 40  # Shoulder to Elbow
L2 = 28  # Elbow to Wrist
L3 = 15  # Wrist to Hand
L4 = 10  # Hand (Finger)
d = 3.5  # Clearance at each joint

# To get the fingertip as close as possible to the shoulder,
# we need a folding configuration that minimizes the final x-coordinate.
# After testing combinations, the (+L1, -L2, -L3, +L4) pattern is optimal.
# The folds are all in the same direction to prevent self-collision.
final_x = L1 - L2 - L3 + L4

# Each of the three joints (Elbow, Wrist, Hand) adds a vertical separation 'd'.
# As the segments stack up, the final vertical position is 3 times the clearance.
final_y = 3 * d

# Calculate the final distance using the Pythagorean theorem
distance = math.sqrt(final_x**2 + final_y**2)

# Output the equation with the numbers plugged in
print(f"The final distance is calculated using the equation:")
print(f"distance = sqrt( (L1 - L2 - L3 + L4)^2 + (3 * d)^2 )")
print(f"distance = sqrt( ({L1} - {L2} - {L3} + {L4})^2 + (3 * {d})^2 )")
print(f"distance = sqrt( ({final_x})^2 + ({final_y})^2 )")
print(f"distance = sqrt( {final_x**2} + {final_y**2} )")
print(f"distance = sqrt( {final_x**2 + final_y**2} )")
print(f"Final calculated distance: ~{distance:.2f} cm")

# The calculated distance 12.62 cm is closest to answer choice A.