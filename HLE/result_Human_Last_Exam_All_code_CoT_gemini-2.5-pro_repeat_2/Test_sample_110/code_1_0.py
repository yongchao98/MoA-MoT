import math

# Define the segment lengths
L1_shoulder_elbow = 40
L2_elbow_wrist = 28
L3_wrist_hand = 15
L4_hand_finger = 10

# Define the base coordinates
S_shoulder = (0, 0)
E_elbow = (L1_shoulder_elbow, 0)

# To satisfy the non-adjacent constraint dist(L1, L3) >= 1 cm,
# we model the arm in a "Z-fold" configuration where L3 is a horizontal
# line at a height of 1 cm.
h = 1.0

# --- Step 1: Calculate the position of the Wrist (W) ---
# The wrist W is at height h=1 and at a distance L2 from the elbow E.
# (W_x - E_x)^2 + (W_y - E_y)^2 = L2^2
# (W_x - 40)^2 + (1 - 0)^2 = 28^2
# We solve for W_x, choosing the smaller value as the arm folds back.
W_x = E_elbow[0] - math.sqrt(L2_elbow_wrist**2 - h**2)
W_y = h
W_wrist = (W_x, W_y)

# --- Step 2: Calculate the position of the Hand (H) ---
# L3 is a horizontal segment of length 15, extending left from the wrist W.
H_x = W_wrist[0] - L3_wrist_hand
H_y = h
H_hand = (H_x, H_y)

# --- Step 3: Calculate the position of the Finger tip (F) ---
# In a Z-fold, L4 is parallel to L2. We find the unit vector for L2
# (from E to W) and apply it to L4 starting from H.
vec_EW = (W_wrist[0] - E_elbow[0], W_wrist[1] - E_elbow[1])
mag_EW = L2_elbow_wrist # The magnitude is the segment length
unit_vec_EW = (vec_EW[0] / mag_EW, vec_EW[1] / mag_EW)

# The vector for L4 (from H to F) is L4 long and in the same direction.
vec_HF = (L4_hand_finger * unit_vec_EW[0], L4_hand_finger * unit_vec_EW[1])
F_x = H_hand[0] + vec_HF[0]
F_y = H_hand[1] + vec_HF[1]
F_finger = (F_x, F_y)

# --- Step 4: Calculate the final distance from Shoulder to Finger ---
# The distance is the magnitude of the vector from S(0,0) to F.
final_distance = math.sqrt(F_finger[0]**2 + F_finger[1]**2)

print("The final position of the finger is calculated based on a Z-fold configuration.")
print(f"The final x-coordinate of the finger is {F_finger[0]:.4f} cm.")
print(f"The final y-coordinate of the finger is {F_finger[1]:.4f} cm.")
print("\nThe final distance is the square root of the sum of the squares of the coordinates:")
print(f"Distance = sqrt( ({F_finger[0]:.4f})^2 + ({F_finger[1]:.4f})^2 )")
print(f"Distance = sqrt( {F_finger[0]**2:.4f} + {F_finger[1]**2:.4f} )")
print(f"Distance = sqrt( {F_finger[0]**2 + F_finger[1]**2:.4f} )")
print(f"Final calculated distance: {final_distance:.4f} cm")
print("\nThis value is closest to answer choice D.")
