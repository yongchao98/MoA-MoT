import math

def get_heesch_numbers():
    """
    Determines the Heesch number for the polygons in the image.

    The Heesch number of a shape P is the maximum number of layers of copies of P
    that can be placed around it without overlapping. If a shape can tile the plane,
    its Heesch number is considered to be infinity.

    1.  Analysis of the Polygon's Shape:
        The polygons in the image are identical. Each is a non-convex 12-sided figure.
        It can be seen as a hexagon modified on a triangular grid. Three of its
        non-adjacent sides have a triangular "bump" added, while the other three
        opposite sides have a triangular "dent" removed.

    2.  Symmetry Property:
        The polygon possesses central symmetry. This means it is unchanged by a
        180-degree rotation about its center. Each bump on one side corresponds
        to a dent on the opposite side. A 180-degree rotation maps these features
        onto each other.

    3.  Tiling Capability:
        A fundamental theorem in geometry states that any centrally symmetric polygon
        can tile the plane (tessellate). The protrusions (bumps) on one tile fit
        perfectly into the indentations (dents) of its neighbors.

    4.  Heesch Number Conclusion:
        Since the polygon can tile the plane, its Heesch number is infinite.
        Because all three polygons are identical, this applies to each of them.
        The problem statement uses the symbol ∞, which we represent as 'infinity'.
    """
    # The shape tiles the plane, so its Heesch number is infinite.
    # The symbol '∞' can be represented as float('inf') in Python, but for clarity
    # and following the conceptual nature of the problem, we will use the string "infinity".
    heesch_number = "infinity"

    # There are three identical polygons.
    number_of_polygons = 3
    results = [heesch_number] * number_of_polygons

    # Print the results separated by commas.
    print(", ".join(results))

get_heesch_numbers()