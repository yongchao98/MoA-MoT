def solve_horse_area():
    """
    Calculates the reachable area for the horse.
    """
    # Step 1: Define constants
    rope_length = 7 / 2
    house_area = 3

    # The vertices of the house polygon
    house_vertices = [
        (-2, 0), (0, 0), (0, -2), (-1, -2), (-1, -1), (-2, -1)
    ]

    # Verify that the entire house is within the rope's reach by checking its vertices
    # A point (x,y) is reachable if |x| + |y| <= rope_length
    for x, y in house_vertices:
        if abs(x) + abs(y) > rope_length:
            # This case does not happen for this problem, but it's a good check.
            print(f"Error: A vertex of the house { (x,y) } is outside the rope's reach.")
            return

    # Step 2: Calculate the area of the total reachable region (a taxi-cab circle, which is a diamond shape)
    # The area of a taxi-cab circle with radius r is 2 * r^2
    total_area = 2 * rope_length**2

    # Step 3: Subtract the area of the house
    # The problem states the horse can't enter the house. Since the geodesic taxi-cab distance
    # from the origin is not affected by this particular obstacle shape, the reachable area is
    # the total area minus the house area.
    final_area = total_area - house_area

    # Step 4: Print the final equation and result
    # The final code needs to output each number in the final equation.
    print(f"The total area the horse could reach without the house is calculated as 2 * (rope length)^2.")
    print(f"The area of the house is given as 3 unit squares.")
    print(f"The reachable area is the total area minus the house area.")
    print(f"Final Calculation:")
    print(f"2 * ({rope_length})^2 - {house_area} = {total_area} - {house_area} = {final_area}")

solve_horse_area()