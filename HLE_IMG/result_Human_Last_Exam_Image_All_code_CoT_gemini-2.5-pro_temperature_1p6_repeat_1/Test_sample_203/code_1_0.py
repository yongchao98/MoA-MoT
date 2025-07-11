def solve_heesch_numbers():
    """
    Determines the Heesch number for each of the three polygons in the image.

    Analysis:
    1.  The top polygon has three pairs of opposite sides that are complementary
        (one side has a "point" and the opposite has a "notch"). This allows it
        to tile the plane. Therefore, its Heesch number is infinity.

    2.  The middle polygon has three adjacent sides with points and three
        adjacent sides with notches. This shape is a known example in tiling theory
        that can be surrounded by exactly one layer of copies of itself, but no more.
        Therefore, its Heesch number is 1.

    3.  The bottom polygon has a more complex arrangement of points and notches.
        It is another known example from the study of tilings (one of Fontaine's tiles)
        that cannot tile the plane but can be surrounded by three layers of its copies.
        Therefore, its Heesch number is 3.
    """

    # Heesch number for the top polygon
    heesch_1 = "infinity"

    # Heesch number for the middle polygon
    heesch_2 = 1

    # Heesch number for the bottom polygon
    heesch_3 = 3

    # Print the answers in order, separated by commas
    print(f"{heesch_1}, {heesch_2}, {heesch_3}")

solve_heesch_numbers()