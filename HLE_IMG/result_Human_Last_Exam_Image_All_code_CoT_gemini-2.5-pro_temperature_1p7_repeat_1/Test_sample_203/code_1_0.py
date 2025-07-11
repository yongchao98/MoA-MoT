def solve_heesch_numbers():
    """
    Determines the Heesch number for the three polygons in the image.

    The Heesch number of a shape is the maximum number of layers of identical
    shapes that can surround it. If a shape can tile the plane, its Heesch
    number is defined as infinity.

    By observing the image, we can see that each of the three purple polygons
    is surrounded by black "holes". The shape of these holes is identical to the
    polygon itself. This demonstrates that each polygon can perfectly tile the plane
    by fitting into the spaces left by its neighbors.

    Therefore, the Heesch number for all three polygons is infinity.
    """

    # We use the string "∞" to represent infinity for the output.
    infinity_symbol = "∞"

    # Heesch number for the top polygon
    heesch_p1 = infinity_symbol

    # Heesch number for the middle polygon
    heesch_p2 = infinity_symbol

    # Heesch number for the bottom polygon
    heesch_p3 = infinity_symbol

    # Print the final answer as a comma-separated list
    print(f"{heesch_p1},{heesch_p2},{heesch_p3}")

solve_heesch_numbers()