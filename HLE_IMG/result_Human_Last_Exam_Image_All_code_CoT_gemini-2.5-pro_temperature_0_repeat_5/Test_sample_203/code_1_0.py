import math

def solve_heesch_number():
    """
    This function determines the Heesch number for the polygons in the image.

    1.  Analysis of the Polygon: The image shows three identical polygons on a triangular grid.
        The shape is a modified hexagon. It has six sides, and they alternate between
        a triangular point jutting out and a triangular notch going in.

    2.  Tiling Property: The image demonstrates that these polygons can tile the plane.
        The pointed part of one polygon fits perfectly into the notched part of an adjacent
        polygon. This pattern can be repeated indefinitely to cover the entire 2D plane
        without any gaps or overlaps.

    3.  Heesch Number Definition: The Heesch number of a shape is the maximum number of
        layers of copies that can surround it. For a shape that tiles the plane, it can be
        surrounded by an infinite number of layers. Therefore, its Heesch number is defined
        as infinity.

    4.  Conclusion: Since the polygon tiles the plane, its Heesch number is infinity.
        As all three polygons in the figure are identical, they all have the same
        Heesch number.
    """
    # In Python, float('inf') represents infinity. We will use the unicode symbol for display.
    heesch_number = "\u221e"  # Unicode for the infinity symbol ∞

    # The problem asks for the Heesch number for each of the three polygons.
    polygon1_heesch = heesch_number
    polygon2_heesch = heesch_number
    polygon3_heesch = heesch_number

    # Print the results separated by commas.
    print(f"{polygon1_heesch}, {polygon2_heesch}, {polygon3_heesch}")

solve_heesch_number()