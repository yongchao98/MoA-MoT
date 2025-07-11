def solve_tiling_periods():
    """
    Solves for the tiling period of four different shapes on a hexagonal grid.
    The period is the minimum number of tiles needed to form a unit cell that
    can tile the plane by translation.
    """

    # Case 1: Points 13, 31, 23
    # These points are the centers of three adjacent hexagons, forming an equilateral triangle.
    # An equilateral triangle needs a reflected copy to form a rhombus, which can tile by translation.
    # The translational unit contains 2 tiles.
    period_1 = 2

    # Case 2: Points 10, 4, 23, 31
    # These points form a general quadrilateral. Any quadrilateral can tile the plane when
    # combined with a 180-degree rotated copy.
    # The translational unit contains 2 tiles.
    period_2 = 2

    # Case 3: Points 5, 15, 17, 19, 21, 7
    # These points are the vertices of a regular hexagon.
    # A regular hexagon can tile the plane by translation alone.
    # The translational unit contains 1 tile.
    period_3 = 1

    # Case 4: A complex 13-sided polygon (4, 5, ..., 13)
    # This shape is complex and asymmetric. It cannot tile by translation alone, nor does it
    # have 3-fold or 6-fold symmetry for higher-order periodic tilings. The most general
    # way for such a tile to form a translational unit is with a transformed copy.
    # The translational unit contains 2 tiles.
    period_4 = 2

    # The problem asks to output each number in the final equation.
    # We interpret this as printing the final result string.
    print(f"{period_1},{period_2},{period_3},{period_4}")

if __name__ == "__main__":
    solve_tiling_periods()