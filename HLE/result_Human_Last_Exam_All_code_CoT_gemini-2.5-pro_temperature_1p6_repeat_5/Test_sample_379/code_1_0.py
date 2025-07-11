import math

# The problem asks us to find the length of a locus of points C on the surface of a cube
# and then compute a ratio. Let r be the side length of the cube.

# Step 1: Find the length of the locus C.
# Through geometric analysis by unfolding the cube, we find that the locus C is composed of
# six separate circular arcs on the faces of the cube.
# The radius for the geodesic distance is given as r.
# Each of the six arcs subtends an angle of pi/3 radians.
num_arcs = 6
angle_per_arc = math.pi / 3

# The total angle for the complete locus is the sum of the angles of these six arcs.
total_angle = num_arcs * angle_per_arc

# The length of a circular arc is radius * angle. The radius of our circle in the unfolded plane is r.
# So, the total length of the locus C is L = r * total_angle.
# L = r * (6 * pi / 3) = r * (2 * pi)
# For the purpose of calculation, we can let r=1, as it will cancel out in the final ratio.
r = 1.0
length_C = r * total_angle


# Step 2: Divide the length of C by 2*pi*r as requested by the problem.
val_2 = 2
val_pi = math.pi
denominator = val_2 * val_pi * r
ratio = length_C / denominator

# Step 3: Convert the final ratio to a whole number percentage.
percentage = int(round(ratio * 100))

# --- Output ---
print("The analysis shows that the total length of the locus C is a sum of 6 arcs.")
print(f"Number of arcs: {num_arcs}")
print(f"Angle per arc: pi/3 radians")
print(f"Total angle: {num_arcs} * (pi/3) = 2*pi radians")
print(f"Total length of C = radius * total_angle = r * (2*pi).")

print("\nThe problem asks for the ratio (Length of C) / (2*pi*r) as a percentage.")
print("\nFinal Equation Breakdown:")
print(f"Numerator (Length of C)   = {val_2} * pi * r")
print(f"Denominator             = {val_2} * pi * r")
# The terms cancel out.
print(f"Ratio = ({val_2}*pi*r) / ({val_2}*pi*r) = {ratio:.1f}")

print("\nConverting the ratio to a whole number percentage:")
print(f"Percentage = {ratio:.1f} * 100%")
print(f"           = {percentage}%")
<<<100>>>