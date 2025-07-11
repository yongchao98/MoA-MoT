def solve_tiling_periods():
    """
    This function determines and prints the periods for the four tiling problems.
    The period is the number of prototiles in a primitive cell that can tile the plane by translation.
    """

    # Period for shape 1: A triangle formed by three hex centers.
    # It requires a 180-degree rotated copy to form a rhombus, which is the primitive cell.
    # The primitive cell contains 2 tiles.
    period_1 = 2

    # Period for shape 2: A rectangle.
    # A rectangle tiles the plane by translation. The primitive cell is the tile itself.
    # The primitive cell contains 1 tile.
    period_2 = 1

    # Period for shape 3: A Type 1 hexagon (one pair of opposite sides are parallel and equal).
    # It tiles the plane by translation. The primitive cell is the tile itself.
    # The primitive cell contains 1 tile.
    period_3 = 1

    # Period for shape 4: A complex polygon derived from three hexagons.
    # This shape, like a trimer of hexagons, does not tile by translation alone.
    # It forms a translational unit cell when paired with a rotated copy.
    # The primitive cell contains 2 tiles.
    period_4 = 2

    # Combine the results into a comma-separated string
    final_answer = f"{period_1},{period_2},{period_3},{period_4}"

    # Print the explanation and the final answer
    print("The periods for the four tiling shapes are calculated as follows:")
    print(f"1) The shape (13, 31, 23) is a triangle. Period = {period_1}")
    print(f"2) The shape (10, 4, 23, 31) is a rectangle. Period = {period_2}")
    print(f"3) The shape (5, 15, 17, 19, 21, 7) is a tileable hexagon. Period = {period_3}")
    print(f"4) The complex shape (4, 5, ..., 13) requires a partner to tile. Period = {period_4}")
    print("\nFinal list of periods:")
    print(final_answer)

solve_tiling_periods()