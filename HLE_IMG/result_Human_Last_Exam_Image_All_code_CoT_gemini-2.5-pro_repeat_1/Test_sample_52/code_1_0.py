# The problem asks for the period of tiling for four shapes on a hexagonal grid.
# The period is the area of the fundamental repeating unit in terms of single hexagons.

# Case 1: A triangle formed by three hexagon centers.
# Area of the triangle is 1/2 hexagon.
# Two triangles form a rhombus of area 1 hexagon, which tiles the plane.
period_1 = 1

# Case 2: A rectangle with an area equal to one hexagon.
# It tiles the plane by translation.
period_2 = 1

# Case 3: A complex hexagon with an area of 3/4 of a hexagon.
# Four such shapes can be combined to form a tile of area 4 * (3/4) = 3.
# This 3-hexagon unit can tile the plane.
period_3 = 3

# Case 4: The points outline the union of the three main hexagons.
# This 3-hexagon shape tiles the plane by translation.
period_4 = 3

# The final answer requires printing the four periods separated by commas.
print(f"{period_1},{period_2},{period_3},{period_4}")