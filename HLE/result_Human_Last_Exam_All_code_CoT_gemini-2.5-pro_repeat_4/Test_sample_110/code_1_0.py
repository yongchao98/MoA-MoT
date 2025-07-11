import math

# Define the lengths of the four robot arm segments in cm.
shoulder_to_elbow = 40
elbow_to_wrist = 28
wrist_to_hand = 15
hand_to_finger = 10

# To find the minimum possible distance between the finger and the shoulder,
# we model the arm in its most folded state, which forms an "S" shape.
# In an ideal scenario, this means the segments are all aligned on a single line,
# with their directions alternating.
# The distance is calculated by starting with the first segment's length,
# subtracting the second, adding the third, and subtracting the fourth.

distance = shoulder_to_elbow - elbow_to_wrist + wrist_to_hand - hand_to_finger

# Print the equation with all the numbers used in the calculation.
print("To find the closest distance, we calculate the position for an ideal S-shaped fold:")
print(f"Final Distance = {shoulder_to_elbow} cm - {elbow_to_wrist} cm + {wrist_to_hand} cm - {hand_to_finger} cm")
print(f"Calculated Ideal Distance = {distance} cm")
print("\nThis ideal distance is a strong estimate. The actual minimum distance will be slightly different due to the physical constraints on the joints, which prevent perfect 180-degree folds. This result is closest to option F.")
