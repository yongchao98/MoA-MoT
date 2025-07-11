import math

# Define the lengths of the robot arm segments in cm
shoulder_to_elbow_L1 = 40
elbow_to_wrist_L2 = 28
wrist_to_hand_L3 = 15
hand_to_finger_L4 = 10

# Define the clearance constraint
# The problem states "3.5 cm clearance on either side of each joint".
# This is interpreted as a total lateral displacement of 2 * 3.5 = 7 cm for each fold.
clearance_per_fold = 2 * 3.5

# To get the finger as close to the shoulder as possible, the arm must fold back on itself.
# We model this as a 2D problem, as a 3D fold would increase the distance.
# The final x-coordinate is determined by the sequence of segment lengths.
# L1 extends in the positive direction. L2 and L3 fold back (negative direction).
# L4 folds forward again to bring the tip closer to the origin.
final_x = shoulder_to_elbow_L1 - elbow_to_wrist_L2 - wrist_to_hand_L3 + hand_to_finger_L4

# The final y-coordinate is the cumulative lateral displacement from the 3 folding joints.
num_folded_joints = 3
final_y = num_folded_joints * clearance_per_fold

# Calculate the final distance from the shoulder (origin) using the Pythagorean theorem.
distance = math.sqrt(final_x**2 + final_y**2)

# Print the final equation with all the numbers
print("Based on the geometric model, the calculation for the distance is:")
print(f"Final X-coordinate = {shoulder_to_elbow_L1} cm - {elbow_to_wrist_L2} cm - {wrist_to_hand_L3} cm + {hand_to_finger_L4} cm = {final_x} cm")
print(f"Final Y-coordinate = {num_folded_joints} folds * {clearance_per_fold} cm/fold = {final_y} cm")
print(f"Distance = sqrt({final_x}^2 + {final_y}^2)")
print(f"Distance = sqrt({final_x**2} + {final_y**2})")
print(f"Distance = sqrt({final_x**2 + final_y**2})")
print(f"Final Distance ≈ {distance:.2f} cm")

# The calculated value of ~22.14 cm is closest to answer choice B.
print("\nThe calculated distance is approximately 22.14 cm, which is closest to answer choice B (~21.82 cm).")

<<<B>>>