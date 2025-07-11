import math

# Define the segment lengths and minimum distance constraint
L1 = 40.0  # Shoulder to Elbow
L2 = 28.0  # Elbow to Wrist
L3 = 15.0  # Wrist to Hand
L4 = 10.0  # Hand to Finger
D_MIN = 3.5 # The minimum distance constraint to use

# Step 1: Calculate the required separation height 'h' for the parallel model
# The separation 'h' must be at least D_MIN for the S1-S3 clearance.
h_for_s1_s3 = D_MIN

# The separation 'h' is also constrained by the S2-S4 clearance.
# In a parallel model, dist(J4, line_S2) = (L3 * h) / L2.
# This distance must be >= D_MIN, so h >= (L2 * D_MIN) / L3.
h_for_s2_s4 = (L2 * D_MIN) / L3

# The actual 'h' must satisfy both constraints.
h = max(h_for_s1_s3, h_for_s2_s4)

# Step 2: Calculate the coordinates of the hand joint (J4)
# y-coordinate of J4 is simply h.
y4 = h
# x-coordinate of J4 is found by chaining the segments.
# x_j3 = L1 - sqrt(L2^2 - h^2)
# x_j4 = x_j3 - L3
x4 = L1 - math.sqrt(L2**2 - h**2) - L3

# Step 3: Calculate the coordinates of the finger tip (P)
# The finger tip P must be at least D_MIN away from S1 (x-axis).
# To minimize distance to the origin, we place it on this boundary.
yp = D_MIN
# The finger tip P is on a circle of radius L4 around J4.
# (xp - x4)^2 + (yp - y4)^2 = L4^2
# Solve for xp: xp = x4 +/- sqrt(L4^2 - (yp - y4)^2)
# We choose the negative root to bring the point closer to the origin (as x4 is negative).
xp = x4 - math.sqrt(L4**2 - (yp - y4)**2)

# Step 4: Calculate the final distance from the shoulder (origin) to the finger tip
final_distance = math.sqrt(xp**2 + yp**2)

# Output the step-by-step calculation
print("--- Model Parameters ---")
print(f"L1 (Shoulder-Elbow) = {L1} cm")
print(f"L2 (Elbow-Wrist)    = {L2} cm")
print(f"L3 (Wrist-Hand)     = {L3} cm")
print(f"L4 (Hand-Finger)    = {L4} cm")
print(f"d_min (Min Dist)    = {D_MIN} cm")
print("\n--- Calculation Steps ---")
print(f"1. Required vertical separation 'h' for S1 || S3 model:")
print(f"   h must be >= {D_MIN}")
print(f"   h must be >= (L2 * d_min) / L3 = ({L2} * {D_MIN}) / {L3} = {h_for_s2_s4:.4f}")
print(f"   Taking the maximum of these: h = {h:.4f} cm")
print("\n2. Position of hand joint J4(x4, y4):")
print(f"   y4 = h = {y4:.4f} cm")
print(f"   x4 = L1 - sqrt(L2^2 - h^2) - L3 = {L1} - sqrt({L2}^2 - {h:.4f}^2) - {L3} = {x4:.4f} cm")
print("\n3. Position of finger tip P(xp, yp):")
print(f"   To minimize distance, set yp on boundary: yp = d_min = {yp:.4f} cm")
print(f"   xp = x4 - sqrt(L4^2 - (yp - y4)^2) = {x4:.4f} - sqrt({L4}^2 - ({yp:.4f} - {y4:.4f})^2) = {xp:.4f} cm")
print("\n--- Final Distance Calculation ---")
print(f"Distance = sqrt(xp^2 + yp^2)")
print(f"Distance = sqrt(({xp:.4f})^2 + ({yp:.4f})^2)")
print(f"Distance = {final_distance:.4f} cm")

# Final check against options
options = {'A': 11.6, 'B': 56.21, 'C': 25.65, 'D': 12.48, 'E': 39.85, 'F': 16.76, 'G': 45.86}
closest_option = min(options.items(), key=lambda item: abs(item[1] - final_distance))
print(f"\nClosest option is {closest_option[0]}: ~{closest_option[1]}")
print(f"<<<D>>>")