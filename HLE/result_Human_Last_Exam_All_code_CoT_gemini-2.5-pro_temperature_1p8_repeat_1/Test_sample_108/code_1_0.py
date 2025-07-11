import math

# Define the robot arm's parameters
L1 = 40.0  # Shoulder to elbow segment length (cm)
L2 = 28.0  # Elbow to wrist segment length (cm)
L3 = 15.0  # Wrist to hand segment length (cm)
L4 = 10.0  # Hand (finger) segment length (cm)

# Define the non-adjacent collision clearance
d_clear = 1.0  # Minimum distance between non-adjacent segments (cm)

# Step 1: Assume the arm curls back such that a later segment is parallel
# to the first segment (L1) at the minimum clearance distance.
# We will find the optimal position for the hand joint (J3).

# The hand joint J3 must lie on a line parallel to L1 at distance d_clear.
# Let's assume this line is y = d_clear.
y_j3 = d_clear

# Step 2: J3 must be reachable by the L2-L3 chain from the elbow joint J1.
# J1 is at (L1, 0).
# The reachable workspace for J3 is an annulus centered at J1.
# To get J3 as close to the origin as possible, we stretch L2 and L3 fully.
# The distance from J1 to J3 is maximized: R = L2 + L3.
R_j3_from_j1 = L2 + L3

# Step 3: Calculate the x-coordinate of J3.
# J3=(x_j3, y_j3) is on a circle around J1=(L1, 0) with radius R.
# (x_j3 - L1)^2 + (y_j3 - 0)^2 = R^2
# (x_j3 - L1)^2 = R^2 - y_j3^2
# x_j3 = L1 +/- sqrt(R^2 - y_j3^2)
# To get closest to the origin, we want the most negative x, so we take the minus sign.
x_j3 = L1 - math.sqrt(R_j3_from_j1**2 - y_j3**2)

# J3 is now at (x_j3, y_j3).

# Step 4: Position the final segment L4 to bring the fingertip Pf closest to the origin.
# To maintain the 1cm clearance from L1, L4 must not go below y=1.
# The closest Pf can be to the origin is also on the line y=d_clear.
# This means L4 is placed horizontally. We try pointing it right.
x_pf = x_j3 + L4
y_pf = y_j3

# Step 5: Calculate the final distance from the origin (shoulder) to the fingertip.
final_distance = math.sqrt(x_pf**2 + y_pf**2)

# Print the final equation with all numbers
x_j3_val = L1 - math.sqrt((L2 + L3)**2 - d_clear**2)
x_pf_val = x_j3_val + L4

print("The calculation for the final fingertip position (x_f, y_f) is as follows:")
print(f"y_f = {d_clear}")
print(f"x_f = (L1 - sqrt((L2 + L3)^2 - d_clear^2)) + L4")
print(f"x_f = ({L1} - sqrt(({L2} + {L3})^2 - {d_clear}^2)) + {L4}")
print(f"x_f = ({L1} - sqrt({L2+L3}^2 - {d_clear**2})) + {L4}")
print(f"x_f = ({L1} - {math.sqrt((L2 + L3)**2 - d_clear**2)}) + {L4}")
print(f"x_f = {x_j3_val} + {L4}")
print(f"x_f = {x_pf_val}")
print("\nThe final distance is calculated using the Pythagorean theorem:")
print("distance = sqrt(x_f^2 + y_f^2)")
print(f"distance = sqrt({x_pf_val}^2 + {y_pf}^2)")
print(f"distance = {final_distance}")