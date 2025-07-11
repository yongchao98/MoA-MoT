# The plan is to construct two perpendicular diameters of the circle.
# The four points where these diameters intersect the circle are the vertices of the inscribed square.
# Let the given points be O (center) and A (on the circumference).

# 1. L: Draw the first diameter.
#    Draw a straight line through the given points O and A.
#    This line is a diameter. It intersects the circle at A and a new point on the opposite side, let's call it B.
print("L")

# 2. C: Draw the first construction circle.
#    To construct a line perpendicular to diameter AB, we use the standard method of creating a perpendicular bisector.
#    Draw a circle centered at point A, with a radius defined by point B (the full diameter length).
print("C")

# 3. C: Draw the second construction circle.
#    Draw a circle centered at point B, with a radius defined by point A (again, the full diameter length).
#    These two circles will intersect at two new points, let's call them P and Q.
print("C")

# 4. L: Draw the second diameter.
#    Draw a straight line connecting the two new intersection points, P and Q.
#    This line is the perpendicular bisector of the diameter AB. It passes through the center O and is therefore the second diameter.
#    This second diameter intersects the original circle at two new points, C and D.
#    The points A, B, C, and D are the vertices of the inscribed square.
print("L")

# The full sequence of commands is LCCL.
# The final answer format requires printing the final sequence in one block.
print("\nFinal sequence:")
print("LCCL")