# The final sequence is composed of two parts:
# 1. Finding the four vertices of the square using the perpendicular bisector method.
#    - L: Draw the first diameter through the center and the given point.
#    - C: Draw a circle from one end of the diameter to the other.
#    - C: Draw a second circle from the other end of the diameter back to the first.
#    - L: Draw a line through the intersections of the two new circles to create the perpendicular diameter.
# 2. Connecting the four vertices to draw the sides of the square.
#    - L, L, L, L: Four lines to connect the four vertices.

construction_sequence = "LCCLLLLL"

print(construction_sequence)