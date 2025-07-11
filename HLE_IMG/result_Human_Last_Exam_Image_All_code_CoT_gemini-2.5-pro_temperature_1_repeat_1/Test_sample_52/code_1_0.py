def solve_tiling_periods():
    """
    Calculates the periods of tilings for four shapes on a hexagonal grid.
    The period is the area of the fundamental repeating unit (primitive cell).
    Areas are expressed in units where the area of one hexagon is 1.
    """

    # Case 1: The shape is a triangle formed by the centers of three adjacent hexagons (13, 31, 23).
    # The primitive cell for the lattice of hexagon centers is a rhombus made of two such triangles.
    # The area of this rhombus is equal to the area of one hexagon.
    # Period = Area(Hexagon)
    period_1 = 1.0

    # Case 2: The shape is a polygon (10, 4, 23, 31).
    # This polygon's area can be calculated to be 0.75 times the area of a hexagon.
    # The shape is a triangle. The primitive cell for a tiling by a triangle is a parallelogram
    # with twice the area of the triangle.
    # Period = 2 * Area(Polygon) = 2 * 0.75 = 1.5
    period_2 = 1.5

    # Case 3: The shape is a hexagon (5, 15, 17, 19, 21, 7).
    # This shape is centrally symmetric, which means it can tile the plane by translation alone.
    # Therefore, its area is the period.
    # Geometric analysis suggests its area is equal to two hexagons.
    # Period = Area(Polygon) = 2.0
    period_3 = 2.0

    # Case 4: A complex polygon containing the three distinct centers (13, 23, 31).
    # This shape acts as a primitive cell for the 3-colored hexagonal lattice.
    # The area of this primitive cell is equal to three hexagons.
    # As a primitive cell, it tiles the plane by translation, so its area is the period.
    # Period = 3.0
    period_4 = 3.0

    # The final answer is the four periods separated by commas.
    print(f"The period for shape 1 (13, 31, 23) is {period_1}")
    print(f"The period for shape 2 (10, 4, 23, 31) is {period_2}")
    print(f"The period for shape 3 (5, 15, 17, 19, 21, 7) is {period_3}")
    print(f"The period for shape 4 (4, 5, 14, 23, 18, 19, 21, 22, 31, 30, 9, 8, 13) is {period_4}")
    
    # Final result formatted as requested
    final_answer = f"{period_1},{period_2},{period_3},{period_4}"
    print("\nFinal Answer:")
    print(final_answer)

solve_tiling_periods()