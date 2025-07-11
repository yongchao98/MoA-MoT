import math

# Define the lengths of the robot arm segments
l1 = 40  # Shoulder to Elbow
l2 = 28  # Elbow to Wrist
l3 = 15  # Wrist to Hand
l4 = 10  # Hand (finger)

# Define the minimum clearance
clearance = 1

# Calculate the length of the first virtual link (P0 to P2)
# This represents the distance from the shoulder to the wrist if L1 and L2 are folded back.
d1 = l1 - l2

# Calculate the length of the second virtual link (P2 to P4)
# This represents the distance from the wrist to the finger if L3 and L4 are folded back.
# We adjust this length by the clearance value to account for the self-collision constraint.
d2 = l3 - l4 - clearance

# Calculate the final distance using the Pythagorean theorem.
# This models a minimal configuration where the virtual links form a right-angled triangle.
# The distance is one side of the triangle, with d1 as the hypotenuse and d2 as the other side.
final_distance = math.sqrt(d1**2 - d2**2)

# Print the equation and the result
print(f"The lengths of the segments are:")
print(f"L1 (Shoulder to Elbow): {l1} cm")
print(f"L2 (Elbow to Wrist): {l2} cm")
print(f"L3 (Wrist to Hand): {l3} cm")
print(f"L4 (Hand to Finger): {l4} cm")
print(f"Minimum clearance: {clearance} cm\n")
print("The calculation is based on a geometric model of the folded arm:")
print(f"First virtual link length (d1) = L1 - L2 = {l1} - {l2} = {d1} cm")
print(f"Second virtual link length (d2) = L3 - L4 - clearance = {l3} - {l4} - {clearance} = {d2} cm")
print(f"Final Distance = sqrt(d1^2 - d2^2) = sqrt({d1}^2 - {d2}^2)")
print(f"Final Distance = sqrt({d1**2} - {d2**2}) = sqrt({d1**2 - d2**2})")
print(f"Final Distance ≈ {final_distance:.2f} cm")

<<<A>>>