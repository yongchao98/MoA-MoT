# The objective is to find the shortest sequence of commands (L for line, C for circle)
# to construct an inscribed square given a circle's center (O) and a point on its circumference (P1).

# Step 1: Draw a Line (L)
# Draw a line through the center O and the point P1. This line creates a diameter,
# identifying a second point P3 on the circle, opposite to P1.
# Points P1 and P3 are two vertices of the square.
command_1 = "L"

# Step 2: Draw a Circle (C)
# To find the other two vertices, we must construct the perpendicular bisector
# of the diameter P1P3. The first step is to draw a circle centered at P1,
# with its radius defined by the point P3.
command_2 = "C"

# Step 3: Draw a second Circle (C)
# Draw another circle centered at P3, with its radius defined by the point P1.
# This circle will have the same size as the one from Step 2.
command_3 = "C"

# Step 4: Draw a final Line (L)
# The two circles from steps 2 and 3 intersect at two new points. A line drawn
# through these two points is the perpendicular bisector of P1P3 and forms the
# second diameter of the original circle. The endpoints of this new diameter
# are the final two vertices of the square.
command_4 = "L"

# The resulting shortest sequence of commands is LCCL.
final_sequence = command_1 + command_2 + command_3 + command_4
print(final_sequence)