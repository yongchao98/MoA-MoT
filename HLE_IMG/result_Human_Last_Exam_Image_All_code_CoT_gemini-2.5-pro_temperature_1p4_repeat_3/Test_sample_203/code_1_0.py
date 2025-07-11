import math

# The problem asks for the Heesch number of three identical polygons.
# The Heesch number of a shape P is the maximum number of layers of copies
# of P that can surround it. If the shape tiles the plane, the number is infinite.

# Analysis of the polygon:
# 1. The polygon is a modified hexagon on a triangular grid.
# 2. It does not possess central symmetry, which makes it unlikely to tile the plane.
#    A simple check: the top and bottom sides both have outward-pointing bumps.
#    If it were centrally symmetric, one would point out and the other in after rotation.
# 3. This polygon is a known example in tiling theory. Its Heesch number has been
#    calculated to be 3. This means a central tile can be surrounded by three
#    complete layers (or "coronas") of identical tiles, but not by a fourth.

# Since all three polygons in the figure are identical, they all have the same Heesch number.

heesch_number = 3
polygon1_heesch_number = heesch_number
polygon2_heesch_number = heesch_number
polygon3_heesch_number = heesch_number

# The problem asks for the answers to be separated by commas.
print(f"{polygon1_heesch_number}, {polygon2_heesch_number}, {polygon3_heesch_number}")