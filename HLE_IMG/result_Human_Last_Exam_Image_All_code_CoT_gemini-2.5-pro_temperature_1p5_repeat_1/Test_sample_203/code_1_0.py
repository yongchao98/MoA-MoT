import math

def solve_heesch_numbers():
    """
    This function determines the Heesch number for the given polygons.

    The polygon shown is a modified hexagon. A regular hexagon can tile the plane.
    This polygon is modified such that one pair of opposite sides has complementary
    shapes (a 'bump' and a 'dent'). This property allows the shape to continue
    tiling the plane by translation, just like a regular hexagon.

    By definition, if a shape can tile the plane, its Heesch number is considered
    to be infinity.

    Since all three polygons in the figure are identical, they all have the same
    Heesch number. We represent infinity using the unicode character '∞'.
    """

    # The Heesch number for a shape that tiles the plane is infinity.
    heesch_number = "∞"

    # There are three identical polygons.
    number_of_polygons = 3

    # Create a list of the answers.
    answers = [heesch_number] * number_of_polygons

    # Print the answers separated by commas, as requested.
    print(",".join(answers))

solve_heesch_numbers()