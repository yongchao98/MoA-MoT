import math

def solve_heesch_number():
    """
    This function determines and prints the Heesch number for the polygons in the image.

    The Heesch number H(P) of a polygon P is the maximum number of layers of
    copies of P that can surround a central copy of P without overlapping.
    If a polygon can tile the plane, its Heesch number is defined as infinity.

    1.  The three polygons in the figure are identical.
    2.  The polygon shape has central symmetry (it is unchanged by a 180-degree rotation
        about its center).
    3.  A well-known theorem states that all centrally symmetric polygons can tile the plane.
        This is because two copies can form a parallelogram, which tiles by translation.
    4.  Since the polygon can tile the plane, it can be surrounded by an infinite number
        of layers.
    5.  Therefore, its Heesch number is infinity.
    """

    # We use the unicode character for infinity for display purposes.
    heesch_number = "∞"

    # All three polygons are identical, so they have the same Heesch number.
    answers = [heesch_number, heesch_number, heesch_number]

    # Print the answers separated by commas, as requested.
    print(", ".join(answers))

solve_heesch_number()