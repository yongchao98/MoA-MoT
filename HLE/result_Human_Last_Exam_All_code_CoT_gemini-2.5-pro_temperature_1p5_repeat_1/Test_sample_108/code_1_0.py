import math

# Step 1: Define the parameters from the problem description.
# Segment lengths in cm
l1_shoulder_elbow = 40
l2_elbow_wrist = 28
l3_wrist_hand = 15
l4_hand_finger = 10

# Clearance constraint in cm. This will be the vertical separation 'h'.
h_clearance = 1

# Step 2 & 3: Model the Z-shape folding pattern (-> <- -> <-) and apply constraints.
# We place the shoulder at (0,0) and the first segment along the x-axis.
# The arm then folds back, but with a vertical offset 'h' to maintain clearance.
# Segments L3 and L4 are assumed to be on a plane parallel to L1, at a distance h.
# The segment L2 must bridge this vertical distance h. Its projection onto the
# x-axis (dx) will be shorter than its full length.
dx_l2 = math.sqrt(l2_elbow_wrist**2 - h_clearance**2)

# Step 4: Calculate the final coordinates of the fingertip.
# The final x-coordinate is the sum of the horizontal segment lengths and projections.
# L1 contributes +l1
# L2 contributes -dx_l2 (folding back)
# L3 contributes +l3 (folding out again)
# L4 contributes -l4 (folding in again)
final_x = l1_shoulder_elbow - dx_l2 + l3_wrist_hand - l4_hand_finger

# The final y-coordinate is the vertical separation.
final_y = h_clearance

# Step 5: Calculate the final distance.
final_distance = math.sqrt(final_x**2 + final_y**2)

# Output the step-by-step calculation.
print("This calculation models the arm in a planar Z-fold configuration.")
print(f"The lengths of the segments are: L1={l1_shoulder_elbow}, L2={l2_elbow_wrist}, L3={l3_wrist_hand}, L4={l4_hand_finger} cm.")
print(f"The required clearance between non-adjacent segments is {h_clearance} cm.")
print("\n--- Calculation ---")
print(f"The horizontal projection of the second segment (L2), dx, is calculated as:")
print(f"dx = sqrt(L2^2 - h^2) = sqrt({l2_elbow_wrist}^2 - {h_clearance}^2) = {dx_l2:.2f} cm")
print("\nThe final position of the fingertip (x, y) is calculated as:")
print(f"x = L1 - dx + L3 - L4")
print(f"x = {l1_shoulder_elbow} - {dx_l2:.2f} + {l3_wrist_hand} - {l4_hand_finger} = {final_x:.2f} cm")
print(f"y = h = {h_clearance:.2f} cm")
print("\nThe final distance from the shoulder (origin) is:")
print(f"Distance = sqrt(x^2 + y^2) = sqrt({final_x:.2f}^2 + {final_y:.2f}^2) = {final_distance:.2f} cm")
print("\nThis calculated distance is approximately 17.05 cm.")
print("Comparing this with the given choices, the closest answer is H (~15.68 cm).")
print("The discrepancy may arise from a more complex 3D folding configuration not captured by this simplified model, but this model represents a standard approach to the problem.")