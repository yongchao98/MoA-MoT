import math

# Let the side length of the cube be r. We can set r=1 for simplicity, as it will cancel out in the final ratio.
r = 1.0

# --- Calculation of the locus length on each face ---

# P is on the edge between the Top and Front faces.

# 1. Length of the arc on the Top face.
# This arc is defined by a circle of radius r centered at P, clipped by the face's side edges.
# The angle subtended by this arc is pi/3 radians (or 60 degrees).
# Arc Length = radius * angle
L_top = r * math.pi / 3

# 2. Length of the arc on the Front face.
# By symmetry with the Top face, this arc has the same length.
L_front = r * math.pi / 3

# 3. Length of arcs on the Left face.
# The left face is adjacent to both the Top and Front faces.
# One arc is formed by paths crossing from the Top face. This is a quarter-circle of radius r. Length = r * pi / 2.
# A second, distinct arc is formed by paths crossing from the Front face. This is also a quarter-circle. Length = r * pi / 2.
L_left = (r * math.pi / 2) + (r * math.pi / 2)

# 4. Length of arcs on the Right face.
# By symmetry with the Left face, the total length of the two arcs on the Right face is the same.
L_right = L_left

# 5. Length on Bottom and Back faces is 0, as they cannot be reached by a path of length r.

# --- Total Length Calculation ---
# Total length is the sum of the lengths on the four relevant faces.
L_total = L_top + L_front + L_left + L_right

# The problem asks us to divide the length of C by 2*pi*r
divisor = 2 * math.pi * r

# --- Final Ratio and Percentage Calculation ---
ratio = L_total / divisor
percentage = ratio * 100
final_answer = round(percentage)

# --- Output the results step-by-step ---
print(f"Let the side length r = {r}")
print(f"The length of the arc on the Top face is r*pi/3 = {L_top:.4f}")
print(f"The length of the arc on the Front face is r*pi/3 = {L_front:.4f}")
print(f"The total length of the two arcs on the Left face is r*pi = {L_left:.4f}")
print(f"The total length of the two arcs on the Right face is r*pi = {L_right:.4f}")
print("---")
print(f"The total length of the locus C = {L_top:.4f} + {L_front:.4f} + {L_left:.4f} + {L_right:.4f} = {L_total:.4f}")
print(f"The value to divide by is 2*pi*r = {divisor:.4f}")
print("---")
print(f"The ratio is L / (2*pi*r) = {L_total:.4f} / {divisor:.4f} = {ratio:.4f}")
print(f"The percentage is the ratio * 100 = {percentage:.2f}%")
print(f"The final answer as a whole number percentage is {final_answer}")

<<<133>>>