# The objective is to construct an inscribed square.
# The vertices of an inscribed square are the endpoints of two perpendicular diameters.

# Let O be the center of the given circle and A be a point on its circumference.

# Step 1: Draw the first diameter.
# (L) Draw a Line through the given points O and A. This line intersects
# the circle at a new point, B. The segment AB is a diameter.
op_1 = "L"

# Step 2: Begin constructing the perpendicular bisector of the diameter AB.
# (C) Draw a Circle with its center at A and a radius equal to the
# length of the diameter AB. Point B is used to define this radius.
op_2 = "C"

# Step 3: Continue the perpendicular bisector construction.
# (C) Draw a Circle with its center at B and a radius equal to the
# length of the diameter BA. Point A is used to define this radius.
op_3 = "C"

# Step 4: Complete the construction.
# The two circles from steps 2 and 3 intersect at two new points, P and Q.
# (L) Draw a Line connecting P and Q. This line is the perpendicular
# bisector of AB, passing through O, and is thus the second diameter.
# The four endpoints of the two diameters are the vertices of the square.
op_4 = "L"

# The shortest sequence of commands is LCCL.
# Let's form the final string representing this sequence.
# Each letter represents a step in our construction process.
final_sequence = op_1 + op_2 + op_3 + op_4

print(final_sequence)