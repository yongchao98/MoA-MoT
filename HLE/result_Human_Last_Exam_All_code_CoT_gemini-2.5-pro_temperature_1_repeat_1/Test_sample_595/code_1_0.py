import math

# Step 1: Define triangle properties
side_leg = 18

# Step 2: Optimal placement strategy
# We place the triangle based on a vector (vx, vy) such that vx^2 + vy^2 = 18^2.
# The optimal orientation is found by maximizing the number of grid lines crossed.
# This occurs when we choose vy to be just over an integer, and vx is determined.
# The best integer pair (ix, iy) with ix^2+iy^2 < 324 is (16,8).
# So we set vy to be slightly larger than 8.
vy = 8.000001
vx = math.sqrt(side_leg**2 - vy**2)

# We can conceptually place vertex A at (e, e) where e is a tiny offset.
# The number of squares crossed by a segment is determined by the integer parts
# of its start and end coordinates.
# Let's calculate the floor of the coordinates of the vertices relative to A.
# A is at (0,0) for this calculation.
# B is at (vx, vy)
# C is at (-vy, vx)

# Step 3: Calculate the floor of the coordinates
f_vx = math.floor(vx)
f_vy = math.floor(vy)
f_Cx = math.floor(-vy) # floor of a negative number
f_Cy = math.floor(vx)

print(f"To maximize the number of squares crossed, we orient the triangle with legs based on a vector (vx, vy).")
print(f"We find that setting vy to be slightly greater than 8 is optimal.")
print(f"Let vy = {vy:.6f}")
print(f"Then vx = sqrt(18^2 - vy^2) = {vx:.6f}")
print(f"The integer parts of the vertex coordinates (relative to vertex A) are:")
print(f"B: (floor(vx), floor(vy)) = ({f_vx}, {f_vy})")
print(f"C: (floor(-vy), floor(vx)) = ({f_Cx}, {f_Cy})")

# Step 4: Calculate the number of squares for each segment
# Number of squares = 1 + |delta_floor_x| + |delta_floor_y|

# Leg AB: from A(0,0) to B(vx, vy)
N_AB = 1 + abs(f_vx - 0) + abs(f_vy - 0)
print("\nCalculating the number of squares crossed by each segment:")
print(f"Leg AB crosses {N_AB} squares. (Equation: 1 + |{f_vx} - 0| + |{f_vy} - 0| = {N_AB})")

# Leg AC: from A(0,0) to C(-vy, vx)
N_AC = 1 + abs(f_Cx - 0) + abs(f_Cy - 0)
print(f"Leg AC crosses {N_AC} squares. (Equation: 1 + |{f_Cx} - 0| + |{f_Cy} - 0| = {N_AC})")

# Hypotenuse BC: from B(vx, vy) to C(-vy, vx)
N_BC = 1 + abs(f_Cx - f_vx) + abs(f_Cy - f_vy)
print(f"Hypotenuse BC crosses {N_BC} squares. (Equation: 1 + |{f_Cx} - {f_vx}| + |{f_Cy} - {f_vy}| = {N_BC})")

# Step 5: Sum the squares and subtract overlaps
total_counted_squares = N_AB + N_AC + N_BC
print(f"\nThe total number of squares counted is {N_AB} + {N_AC} + {N_BC} = {total_counted_squares}.")

# The squares containing the three vertices A, B, C are each counted twice.
# These 3 squares are the only overlaps for this generic orientation.
num_overlaps = 3
print(f"The squares containing the three vertices are counted in two segment calculations each. So we subtract {num_overlaps} overlaps.")

k = total_counted_squares - num_overlaps
print(f"\nThe largest number of grid squares, k, is {total_counted_squares} - {num_overlaps} = {k}.")
print(f"The final equation is: {k} = {N_AB} + {N_AC} + {N_BC} - {num_overlaps}")