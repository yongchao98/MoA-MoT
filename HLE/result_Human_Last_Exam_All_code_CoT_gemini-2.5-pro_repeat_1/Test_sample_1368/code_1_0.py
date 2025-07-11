# This script provides the solution to the geometric construction puzzle.
# The goal is to find the shortest sequence of commands (L for Line, C for Circle)
# to construct a square inscribed in a given circle.

# The derivation of the sequence is based on the following geometric principles:
# 1. (L) A line through the center and a given point on the circumference defines the first diameter.
# 2. (C, C) Two circles, centered on the endpoints of the first diameter and with a radius equal to the
#    diameter's length, are drawn.
# 3. (L) A line connecting the intersection points of these two circles creates the second diameter,
#    which is perpendicular to the first.
# This sequence of constructions successfully defines the four vertices of the inscribed square.

shortest_sequence = "LCCL"

print(shortest_sequence)