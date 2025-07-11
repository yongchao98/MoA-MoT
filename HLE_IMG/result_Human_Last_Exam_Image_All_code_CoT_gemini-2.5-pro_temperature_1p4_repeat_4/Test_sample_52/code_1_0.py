def solve():
    """
    Calculates the periods of the four tilings described.
    """
    # Case 1: Triangle of hex centers
    # Area of tile is 1/2 area of a hexagon.
    # Unit cell is a parallelogram made of 2 tiles, so area is 2 * 1/2 = 1 hexagon area.
    period1 = 1

    # Case 2: Quadrilateral 10-4-23-31
    # Area of tile is 9/8 area of a hexagon.
    # To form a repeating unit cell, we need k tiles such that k * 9/8 = M (integer).
    # Smallest k is 8, giving a unit cell area M = 9 hexagon areas.
    period2 = 9

    # Case 3: Hexagon 5-15-17-19-21-7
    # This forms a large hexagon with an area equivalent to 4 small hexagons.
    # This large hexagon tiles the plane, so the unit cell is one tile.
    # The area of the unit cell is 4 hexagon areas.
    period3 = 4

    # Case 4: Complex polygon involving points from hexagons H1, H2, and H3
    # The polygon encloses the three central hexagons.
    # It's assumed its area is the sum of the areas of these three hexagons.
    # Area of the unit cell is 3 hexagon areas.
    period4 = 3
    
    print(f"{period1},{period2},{period3},{period4}")

solve()