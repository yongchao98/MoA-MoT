def solve_heesch_number():
    """
    Solves for the Heesch number of the polygons in the image.

    Step-by-step reasoning:
    1.  The image displays three purple polygons against a triangular grid.
    2.  Careful visual inspection reveals that all three polygons are identical in shape. Therefore, they must have the same Heesch number.
    3.  The key geometric property of this polygon is its symmetry. The shape has 180-degree rotational symmetry (point symmetry or central symmetry) around its center. If rotated by 180 degrees, it maps onto itself.
    4.  A well-established theorem in geometry states that any centrally symmetric polygon tiles the plane. Such polygons are a class of shapes known as parallelogons.
    5.  The problem states that if a polygon tiles the plane, its Heesch number is defined to be infinity. This is because a tiling represents an infinite number of layers surrounding any given tile.
    6.  Since the polygon is centrally symmetric, it tiles the plane, and thus its Heesch number is infinity.
    7.  As all three polygons are the same shape, the answer is "infinity" for each one.
    """

    # The Heesch number for a polygon that tiles the plane is infinity.
    heesch_number = "infinity"

    # Since all three polygons are identical and tile the plane, they all have the same Heesch number.
    results = [heesch_number, heesch_number, heesch_number]

    # The final answer should be the three Heesch numbers separated by commas.
    print(", ".join(results))

solve_heesch_number()