import math

def get_heesch_numbers():
    """
    This function determines the Heesch number for the polygons in the image.

    1.  Analysis of the Polygon: The polygon shown is a modified hexagon on a
        triangular grid. It has alternating protrusions and indentations.
        All three polygons in the figure are identical.

    2.  Tiling Property: A shape's Heesch number is defined as infinity (∞) if
        it can tile the plane. We test if this polygon can tile the plane.
        - The protrusions are shaped to fit perfectly into the indentations.
        - The image shows that the polygons can be stacked vertically by simple
          translation, with the bottom of one fitting the top of the next.
        - The shape can also be translated diagonally (up and to the right) to
          fit against itself, allowing columns of tiles to fill the plane.
        - Since the polygon can tile the plane using two translation vectors,
          it is a "tiler".

    3.  Conclusion: As the polygon tiles the plane, its Heesch number is ∞.
        Since all three polygons are identical, they all have a Heesch number
        of ∞.
    """
    # In mathematics, if a shape tiles the plane, its Heesch number is infinity.
    # Python's float('inf') represents infinity, but we will use the unicode
    # symbol '∞' for the output as requested by the problem's context.
    heesch_number_for_tiler = '∞'

    # The problem shows three identical polygons.
    num_polygons = 3
    
    # The Heesch number is the same for all three.
    heesch_numbers = [heesch_number_for_tiler] * num_polygons
    
    # Print the result as a comma-separated string.
    print(", ".join(heesch_numbers))

get_heesch_numbers()