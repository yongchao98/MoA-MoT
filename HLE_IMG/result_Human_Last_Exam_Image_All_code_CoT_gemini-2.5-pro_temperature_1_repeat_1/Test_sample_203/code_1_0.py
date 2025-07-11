import math

def solve_heesch_number():
    """
    This function determines the Heesch number for the polygons in the image.

    Step 1: Analyze the polygon's geometry.
    The polygon is a dodecagon (12-sided) derived from a hexagon.
    Its boundary is formed by alternating "bumps" (added equilateral triangles)
    and "dents" (removed equilateral triangles) on the sides of a base hexagon.

    Step 2: Determine if the polygon can tile the plane.
    A key property of this polygon is that each pair of its opposite sides has
    complementary modifications. One side has a bump, and the opposite side has
    a matching dent. This allows the shape to tile the plane seamlessly through
    translation, as the bump of one tile fits into the dent of an adjacent tile.

    Step 3: Define the Heesch number for a tiling polygon.
    The Heesch number is the maximum number of layers of a polygon that can
    surround a central one. If a polygon can tile the plane, it can be surrounded
    by an infinite number of layers. By definition, its Heesch number is infinity.

    Step 4: Conclude for all three polygons.
    The image shows three identical polygons. Since the shape tiles the plane,
    each polygon has a Heesch number of infinity.
    """

    # In mathematics, the symbol for infinity is ∞. We represent it as a string.
    # We could also use float('inf') but the symbol is more conventional here.
    heesch_number = "∞"

    # There are three identical polygons in the figure.
    num_polygons = 3

    # The Heesch number is the same for all three.
    results = [heesch_number] * num_polygons

    # Print the answers in order, separated by commas.
    print(','.join(results))

solve_heesch_number()