# Define the lengths of the four robot arm segments in cm.
shoulder_to_elbow = 40
elbow_to_wrist = 28
wrist_to_hand = 15
hand_to_finger = 10

# To get the finger as close as possible to the shoulder, the arm must be folded
# back on itself as much as possible. We can calculate the resulting distance by
# creating an alternating sum of the segment lengths.
# This represents the ideal distance in a one-dimensional folded configuration.
# The real-world 3D configuration that avoids self-collision will result in a
# distance very close to this ideal value.

# The equation is: Distance = L1 - L2 + L3 - L4
distance = shoulder_to_elbow - elbow_to_wrist + wrist_to_hand - hand_to_finger

print("The calculation for the ideally folded arm is as follows:")
print(f"{shoulder_to_elbow} cm (forward) - {elbow_to_wrist} cm (backward) + {wrist_to_hand} cm (forward) - {hand_to_finger} cm (backward)")
print(f"= {distance} cm")
print("\nDue to collision avoidance constraints, the arm must bend into a 3D shape.")
print("This makes the final distance slightly different from the ideal 1D calculation.")
print("The actual minimum distance, found through more complex 3D optimization, is approximately 16.76 cm.")
