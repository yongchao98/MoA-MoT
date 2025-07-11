def solve_periods():
    """
    This function prints the calculated periods for the four tiling problems.
    The reasoning is as follows:
    The period is the number of hexagons in the minimal repeating translational unit cell.

    1) The motif {13, 31, 23} can form a rhombus tile. A valid tiling using this tile has a unit cell containing 2 base hexagons. Period = 2.
    2) The motif {10, 4, 23, 31} is a rectangle. A valid tiling can be constructed where the unit cell contains 2 base hexagons. Period = 2.
    3) The motif {5, 15, 17, 19, 21, 7} is simply a hexagon. Tiling with hexagons has a period of 1.
    4) The motif is a complex shape. To form a tile that can create a valid tiling on the grid, 3 copies of the motif must be combined. This composite tile has an area equivalent to 4 base hexagons. Period = 4.
    """
    period1 = 2
    period2 = 2
    period3 = 1
    period4 = 4

    # Phrase the answer as the four periods separated by ","
    print(f"{period1},{period2},{period3},{period4}")

solve_periods()