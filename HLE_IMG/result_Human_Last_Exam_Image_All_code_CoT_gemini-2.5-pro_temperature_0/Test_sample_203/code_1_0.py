import math

def solve_heesch_number():
    """
    This function determines and prints the Heesch number for the polygons in the image.

    The Heesch number of a shape is the maximum number of layers of congruent
    copies of the shape that can surround it. If a shape can tile the plane,
    its Heesch number is defined as infinity.

    1.  The polygon in the image is a non-convex, centrally symmetric hexagon.
    2.  A mathematical theorem states that all centrally symmetric polygons can tile the plane.
    3.  Therefore, the polygon shown can tile the plane.
    4.  According to the problem's definition, its Heesch number is infinity.
    5.  Since all three polygons in the figure are identical, they all have a Heesch number of infinity.
    """

    # Representing infinity. In Python, float('inf') is the standard representation.
    # For clear textual output as requested, we will use the string "infinity".
    heesch_number_1 = "infinity"
    heesch_number_2 = "infinity"
    heesch_number_3 = "infinity"

    # The final answer is the Heesch number for each of the three polygons, separated by commas.
    print(f"{heesch_number_1}, {heesch_number_2}, {heesch_number_3}")

solve_heesch_number()