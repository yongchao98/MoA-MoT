import math

# Define the segment lengths and clearance
L1 = 40.0  # Shoulder to Elbow
L2 = 28.0  # Elbow to Wrist
L3 = 15.0  # Wrist to Hand
L4 = 10.0  # Hand to Finger
h = 3.5   # Clearance at joints

print(f"Robot arm segment lengths:")
print(f"L1 (Shoulder-Elbow): {L1} cm")
print(f"L2 (Elbow-Wrist): {L2} cm")
print(f"L3 (Wrist-Hand): {L3} cm")
print(f"L4 (Hand-Finger): {L4} cm")
print(f"Joint clearance (h): {h} cm\n")

# Step 1: Position Shoulder and Elbow
S = (0.0, 0.0)
E = (L1, 0.0)
print(f"Positioning Shoulder S at {S}")
print(f"Positioning Elbow E at {E}")

# Step 2: Calculate Wrist (W) position
# We fold L2 back towards L1. This creates a triangle SEW.
# The height from W to the line SE is h=3.5.
# We use the Pythagorean theorem to find the horizontal projection of L2.
x_proj_L2 = math.sqrt(L2**2 - h**2)
# The x-coordinate of W is E.x minus this projection.
x_W = E[0] - x_proj_L2
# The y-coordinate of W is h (we fold "up").
y_W = h
W = (x_W, y_W)
print(f"\nFolding at Elbow joint:")
print(f"The horizontal projection of the Elbow-Wrist segment is sqrt({L2}^2 - {h}^2) = {x_proj_L2:.4f} cm")
print(f"Positioning Wrist W at ({E[0]} - {x_proj_L2:.4f}, {h}) = {W}")


# To calculate the next positions, we use local coordinate frames and rotations.
def rotate_and_translate(point_rel, angle_rad, origin_abs):
    """Rotates a relative point and translates it to an absolute origin."""
    cos_a = math.cos(angle_rad)
    sin_a = math.sin(angle_rad)
    x_abs = point_rel[0] * cos_a - point_rel[1] * sin_a + origin_abs[0]
    y_abs = point_rel[0] * sin_a + point_rel[1] * cos_a + origin_abs[1]
    return (x_abs, y_abs)

# Step 3: Calculate Hand (H) position
# We fold L3 back towards L2. To achieve the tightest curl, we fold "down".
# In the local frame of W (with x-axis along WE), H has coordinates (-x_proj, -h).
x_proj_L3 = math.sqrt(L3**2 - h**2)
H_rel = (-x_proj_L3, -h)
# We need to rotate this relative position back to the main coordinate system.
# The rotation angle is the angle of the vector WE.
vec_WE = (E[0] - W[0], E[1] - W[1])
angle_WE = math.atan2(vec_WE[1], vec_WE[0])
H = rotate_and_translate(H_rel, angle_WE, W)
print(f"\nFolding at Wrist joint:")
print(f"The horizontal projection of the Wrist-Hand segment is sqrt({L3}^2 - {h}^2) = {x_proj_L3:.4f} cm")
print(f"Positioning Hand H at {H}")

# Step 4: Calculate Finger (F) position
# We fold L4 back towards L3. We continue the "downward" curl.
# In the local frame of H (with x-axis along HW), F has coordinates (-x_proj, -h).
x_proj_L4 = math.sqrt(L4**2 - h**2)
F_rel = (-x_proj_L4, -h)
# The rotation angle is the angle of the vector HW.
vec_HW = (W[0] - H[0], W[1] - H[1])
angle_HW = math.atan2(vec_HW[1], vec_HW[0])
F = rotate_and_translate(F_rel, angle_HW, H)
print(f"\nFolding at Hand joint:")
print(f"The horizontal projection of the Hand-Finger segment is sqrt({L4}^2 - {h}^2) = {x_proj_L4:.4f} cm")
print(f"Positioning Finger F at {F}")

# Step 5: Calculate the final distance from Shoulder (S) to Finger (F)
distance = math.sqrt(F[0]**2 + F[1]**2)
print(f"\nFinal Calculation:")
print(f"The final distance between the Shoulder S(0,0) and Finger F{F} is:")
print(f"distance = sqrt({F[0]:.4f}^2 + {F[1]:.4f}^2) = {distance:.4f} cm")
