import math

# Define segment lengths and separation distance
L1 = 40  # Shoulder to Elbow
L2 = 28  # Elbow to Wrist
L3 = 15  # Wrist to Hand
L4 = 10  # Hand (Finger)

# From "3.5 cm clearance on either side", we assume a total horizontal separation
# between the centerlines of adjacent folded segments.
dx = 3.5 + 3.5

# --- Calculations ---

# Initial joint positions
p0_x, p0_y = 0, 0

# P1 (Elbow) - Segment L1 extends straight up
p1_x = p0_x
p1_y = p0_y + L1

# P2 (Wrist) - Segment L2 folds back (dy is negative)
# dy2 = sqrt(L2^2 - dx^2)
dy2 = math.sqrt(L2**2 - dx**2)
p2_x = p1_x + dx
p2_y = p1_y - dy2

# P3 (Hand) - Segment L3 folds back again (dy is negative)
# dy3 = sqrt(L3^2 - dx^2)
dy3 = math.sqrt(L3**2 - dx**2)
p3_x = p2_x + dx
p3_y = p2_y - dy3

# P4 (Finger Tip) - Segment L4 folds forward (dy is positive)
# dy4 = sqrt(L4^2 - dx^2)
dy4 = math.sqrt(L4**2 - dx**2)
p4_x = p3_x + dx
p4_y = p3_y + dy4

# Final distance from origin
final_distance = math.sqrt(p4_x**2 + p4_y**2)

# --- Output ---
print("Step-by-step calculation of the final finger tip position (P4):")
print(f"Segment lengths (cm): L1={L1}, L2={L2}, L3={L3}, L4={L4}")
print(f"Assumed horizontal separation between segments (dx) = {dx:.2f} cm\n")

print("1. Position of Elbow (P1):")
print(f"   P1 = (x: {p0_x}, y: {p0_y}) + (dx: 0, dy: {L1}) = ({p1_x:.2f}, {p1_y:.2f}) cm\n")

print("2. Position of Wrist (P2):")
print(f"   Vertical component of L2 (dy2) = sqrt({L2}^2 - {dx}^2) = {dy2:.2f} cm")
print(f"   P2 = (x: {p1_x:.2f}, y: {p1_y:.2f}) + (dx: {dx}, dy: -{dy2:.2f}) = ({p2_x:.2f}, {p2_y:.2f}) cm\n")

print("3. Position of Hand (P3):")
print(f"   Vertical component of L3 (dy3) = sqrt({L3}^2 - {dx}^2) = {dy3:.2f} cm")
print(f"   P3 = (x: {p2_x:.2f}, y: {p2_y:.2f}) + (dx: {dx}, dy: -{dy3:.2f}) = ({p3_x:.2f}, {p3_y:.2f}) cm\n")

print("4. Position of Finger Tip (P4):")
print(f"   Vertical component of L4 (dy4) = sqrt({L4}^2 - {dx}^2) = {dy4:.2f} cm")
print(f"   P4 = (x: {p3_x:.2f}, y: {p3_y:.2f}) + (dx: {dx}, dy: +{dy4:.2f}) = ({p4_x:.2f}, {p4_y:.2f}) cm\n")

print("Final Distance Calculation:")
print(f"Distance = sqrt(P4_x^2 + P4_y^2)")
print(f"Distance = sqrt({p4_x:.2f}^2 + {p4_y:.2f}^2)")
print(f"Distance = sqrt({p4_x**2:.2f} + {p4_y**2:.2f})")
print(f"Distance = sqrt({(p4_x**2 + p4_y**2):.2f})")
print(f"Final distance = {final_distance:.2f} cm")
