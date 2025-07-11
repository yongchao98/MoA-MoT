def find_tiling_periods():
    """
    This function prints the periods of four different tile patterns on a hexagonal grid.
    The periods are determined based on the geometric properties of the shapes and tiling theory.
    """
    # Period for shape 1) 13, 31, 23 (Asymmetric triangle)
    period_1 = 2

    # Period for shape 2) 10, 4, 23, 31 (Isosceles trapezoid)
    period_2 = 2

    # Period for shape 3) 5, 15, 17, 19, 21, 7 (Type 1 space-filling hexagon)
    period_3 = 1

    # Period for shape 4) 4, 5, ..., 13 (Complex asymmetric polygon)
    period_4 = 2

    # The final answer is the four periods separated by a comma.
    print(f"{period_1},{period_2},{period_3},{period_4}")

find_tiling_periods()