import math

# This problem is solved by understanding the geometry and symmetry of an icosahedron.

# 1. An icosahedron standing on a face has a parallel face at the top.
# 2. When half-filled, the water level is at the horizontal plane exactly
#    midway between the top and bottom faces, as this plane bisects the volume.
# 3. This mid-plane intersects 6 of the icosahedron's edges.
# 4. Due to the icosahedron's symmetry, these 6 intersection points form a regular polygon.

# Let's define the final shape with an equation for its number of sides.
# Let S be the number of sides of the polygon.
S = 6
shape_name = "regular hexagon"

# Output the reasoning and the final answer.
print("The shape of the water surface is a polygon whose properties can be described by the following equation:")
print(f"Number of sides (S) = {S}")
print(f"A polygon with {S} sides is a hexagon.")
print(f"Because of the icosahedron's symmetry, the shape is a {shape_name}.")
