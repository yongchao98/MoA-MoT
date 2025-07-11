def calculate_all_periods():
    """
    Calculates the tiling period for each of the four specified shapes on the hexagonal grid.
    The period is the number of tiles in the smallest translational unit cell.
    """

    # Case 1: An equilateral triangle formed by three hexagon centers.
    # Unit cell is a rhombus containing 2 tiles.
    period_1 = 2

    # Case 2: A rectangle.
    # Unit cell is the rectangle itself, containing 1 tile.
    period_2 = 1

    # Case 3: A half-hexagon.
    # Unit cell is a full hexagon containing 2 tiles.
    period_3 = 2

    # Case 4: A three-hexagon cluster.
    # Unit cell is a six-hexagon block containing 2 tiles.
    period_4 = 2

    # The final answer is the four periods separated by commas.
    final_answer = f"{period_1},{period_2},{period_3},{period_4}"
    print(final_answer)

calculate_all_periods()